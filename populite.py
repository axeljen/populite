#!/usr/bin/env python3
"""
Mask-aware windowed and genome-wide population statistics from a SNP-only
VCF plus per-sample callability BED masks. Select outputs with --stats; the
default is pi and het.

pi (per population) follows the pixy estimator: the sum of pairwise haplotype
differences divided by the number of pairwise haplotype comparisons, where
invariant callable sites are counted from the per-sample masks.

Heterozygosity (per sample) is the number of heterozygous calls inside the
sample's mask divided by the number of callable bases in its mask (minus
variant sites inside the mask where the genotype is missing).

dxy is the mean number of cross-population haplotype differences per
callable cross-population comparison. FST uses Weir and Cockerham's variance
components for diploid, biallelic SNPs. Watterson's theta and Tajima's D use
per-site called-allele counts to account for missing genotypes.
"""

import argparse
import gzip
import itertools
import subprocess
import sys
from collections import defaultdict

import numpy as np


# ----------------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------------

def read_chroms_file(path):
    """Read chromosome names, ignoring blank lines and full-line comments."""
    names = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            if len(name.split()) != 1 or "," in name:
                raise ValueError(
                    f"{path}:{lineno}: expected one chromosome name per row"
                )
            names.append(name)
    return names


def read_two_col(path):
    """Read a two-column whitespace-separated file into an ordered list of pairs."""
    pairs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                sys.exit(f"ERROR: expected two columns in {path}: {line!r}")
            pairs.append((fields[0], fields[1]))
    return pairs


def merge_intervals(starts, ends):
    """Sort and merge overlapping/abutting half-open intervals."""
    if len(starts) == 0:
        return starts, ends
    order = np.argsort(starts, kind="stable")
    starts, ends = starts[order], ends[order]
    cummax = np.maximum.accumulate(ends)
    new = np.empty(len(starts), dtype=bool)
    new[0] = True
    new[1:] = starts[1:] > cummax[:-1]
    idx = np.flatnonzero(new)
    return starts[idx], np.maximum.reduceat(ends, idx)


def read_bed(path):
    """Read a BED file into {chrom: (starts, ends)} with merged intervals."""
    opener = gzip.open if path.endswith(".gz") else open
    raw = defaultdict(lambda: ([], []))
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t", 3)[:3]
            s, e = raw[chrom]
            s.append(int(start))
            e.append(int(end))
    out = {}
    for chrom, (s, e) in raw.items():
        out[chrom] = merge_intervals(np.asarray(s, dtype=np.int64),
                                     np.asarray(e, dtype=np.int64))
    return out


# ----------------------------------------------------------------------------
# Window helpers
# ----------------------------------------------------------------------------

def weighted_window_lengths(starts, ends, weights, wsize, nwin):
    """
    Sum of (overlap length * weight) per window for half-open intervals.
    Intervals may span several windows.
    """
    out = np.zeros(nwin, dtype=np.float64)
    if len(starts) == 0:
        return out
    ws = starts // wsize
    we = (ends - 1) // wsize
    same = ws == we
    np.add.at(out, ws[same], (ends[same] - starts[same]) * weights[same])
    span = ~same
    if span.any():
        s, e, w = starts[span], ends[span], weights[span]
        a, b = ws[span], we[span]
        np.add.at(out, a, ((a + 1) * wsize - s) * w)
        np.add.at(out, b, (e - b * wsize) * w)
        # full windows strictly between a and b
        diff = np.zeros(nwin + 1, dtype=np.float64)
        np.add.at(diff, a + 1, w)
        np.add.at(diff, b, -w)
        out += np.cumsum(diff[:-1]) * wsize
    return out


def clip(starts, ends, chrom_len):
    ends = np.minimum(ends, chrom_len)
    keep = ends > starts
    return starts[keep], ends[keep]


def coverage_segments(interval_sets):
    """
    Sweep-line over several sets of merged intervals. Returns segment
    starts, ends and the number of sets covering each segment (>0 only).
    """
    pos = np.concatenate([np.concatenate([s, e]) for s, e in interval_sets])
    delta = np.concatenate([np.concatenate([np.ones(len(s), np.int64),
                                            -np.ones(len(e), np.int64)])
                            for s, e in interval_sets])
    if len(pos) == 0:
        empty = np.zeros(0, np.int64)
        return empty, empty, empty
    upos, inv = np.unique(pos, return_inverse=True)
    depth = np.cumsum(np.bincount(inv, weights=delta).astype(np.int64))
    seg_s, seg_e, seg_n = upos[:-1], upos[1:], depth[:-1]
    keep = seg_n > 0
    return seg_s[keep], seg_e[keep], seg_n[keep]


def choose2(x):
    x = np.asarray(x, dtype=np.float64)
    return x * (x - 1) / 2


def tajima_constants(n):
    """Return a1, e1 and e2 for Tajima's D at sample size n."""
    if n < 2:
        return None, None, None
    i = np.arange(1, n, dtype=np.float64)
    a1 = np.sum(1 / i)
    a2 = np.sum(1 / (i * i))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    return a1, c1 / a1, c2 / (a1 * a1 + a2)


def tajima_stdev(total_alleles, n_sites, mutations):
    """Pixy-compatible Tajima's D denominator with site-varying missingness."""
    if n_sites <= 0 or mutations <= 0:
        return None
    n = int(np.rint(total_alleles / n_sites))
    _, e1, e2 = tajima_constants(n)
    if e1 is None:
        return None
    variance = e1 * mutations + e2 * mutations * (mutations - 1)
    return np.sqrt(variance) if variance > 0 else None


def wc_fst_components(G, pop1, pop2, max_allele):
    """Return per-site Weir-Cockerham a, b, c sums for diploid genotypes."""
    G1 = G[:, pop1, :]
    G2 = G[:, pop2, :]
    allele_counts = []
    counts1, counts2 = [], []
    for allele in range(max_allele + 1):
        allele_counts.append((G == allele).sum(axis=(1, 2)))
        counts1.append((G1 == allele).sum(axis=(1, 2)))
        counts2.append((G2 == allele).sum(axis=(1, 2)))
    allele_counts = np.column_stack(allele_counts)
    counts1, counts2 = np.column_stack(counts1), np.column_stack(counts2)
    biallelic = (allele_counts > 0).sum(axis=1) == 2

    an1, an2 = counts1.sum(axis=1), counts2.sum(axis=1)
    n1, n2 = an1 // 2, an2 // 2
    ntotal = n1 + n2
    nbar = ntotal / 2
    nc = ntotal - (n1 * n1 + n2 * n2) / np.maximum(ntotal, 1)
    valid = biallelic & (an1 > 0) & (an2 > 0) & (nbar > 1) & (nc > 0)

    a_site = np.zeros(len(G), dtype=np.float64)
    b_site = np.zeros(len(G), dtype=np.float64)
    c_site = np.zeros(len(G), dtype=np.float64)
    if not valid.any():
        return a_site, b_site, c_site, valid

    p1 = counts1 / np.maximum(an1[:, None], 1)
    p2 = counts2 / np.maximum(an2[:, None], 1)
    pbar = (counts1 + counts2) / np.maximum((an1 + an2)[:, None], 1)
    s_squared = (n1[:, None] * (p1 - pbar) ** 2 +
                 n2[:, None] * (p2 - pbar) ** 2) / np.maximum(nbar[:, None], 1)

    het_counts = []
    for allele in range(max_allele + 1):
        h1 = ((G1[:, :, 0] >= 0) & (G1[:, :, 1] >= 0) &
              (G1[:, :, 0] != G1[:, :, 1]) &
              ((G1[:, :, 0] == allele) | (G1[:, :, 1] == allele))).sum(axis=1)
        h2 = ((G2[:, :, 0] >= 0) & (G2[:, :, 1] >= 0) &
              (G2[:, :, 0] != G2[:, :, 1]) &
              ((G2[:, :, 0] == allele) | (G2[:, :, 1] == allele))).sum(axis=1)
        het_counts.append(h1 + h2)
    hbar = np.column_stack(het_counts) / np.maximum(ntotal[:, None], 1)

    for allele in range(max_allele + 1):
        p = pbar[:, allele]
        s2 = s_squared[:, allele]
        h = hbar[:, allele]
        within = p * (1 - p) - s2 / 2
        a = (nbar / np.maximum(nc, 1)) * (
            s2 - (within - h / 4) / np.maximum(nbar - 1, 1))
        b = (nbar / np.maximum(nbar - 1, 1)) * (
            within - ((2 * nbar - 1) * h / (4 * np.maximum(nbar, 1))))
        c = h / 2
        a_site += np.where(valid, a, 0)
        b_site += np.where(valid, b, 0)
        c_site += np.where(valid, c, 0)
    return a_site, b_site, c_site, valid


def intersect_intervals(starts1, ends1, starts2, ends2):
    """Return intersections of two sorted, merged half-open interval sets."""
    starts, ends = [], []
    i = j = 0
    while i < len(starts1) and j < len(starts2):
        start = max(starts1[i], starts2[j])
        end = min(ends1[i], ends2[j])
        if end > start:
            starts.append(start)
            ends.append(end)
        if ends1[i] <= ends2[j]:
            i += 1
        else:
            j += 1
    return np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)


# ----------------------------------------------------------------------------
# Variant reading
# ----------------------------------------------------------------------------

def read_vcf_header(vcf_path):
    """Return (samples, {contig: length or None}) from the VCF header via bcftools."""
    res = subprocess.run(["bcftools", "view", "-h", vcf_path],
                         stdout=subprocess.PIPE, universal_newlines=True)
    if res.returncode != 0:
        sys.exit(f"ERROR: could not read VCF header of {vcf_path}")
    samples, contigs = [], {}
    for line in res.stdout.splitlines():
        if line.startswith("##contig=<"):
            fields = dict(kv.split("=", 1) for kv in line[10:].rstrip(">").split(",") if "=" in kv)
            length = fields.get("length")
            contigs[fields["ID"]] = int(length) if length else None
        elif line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
    return samples, contigs


def parse_gt_chunk(gt_lines, nsamp, ploidy):
    """Parse tab-separated GT strings (one line per site) into int8 (k, nsamp, ploidy)."""
    k = len(gt_lines)
    text = "\t".join(gt_lines).replace("|", "/").replace(".", "-1")
    toks = text.replace("/", "\t").split("\t")
    if len(toks) == k * nsamp * ploidy:
        return np.array(toks).astype(np.int8).reshape(k, nsamp, ploidy)
    # ploidy differs for some genotypes: parse one by one
    G = np.full((k, nsamp, ploidy), -1, dtype=np.int8)
    for i, line in enumerate(gt_lines):
        for j, gt in enumerate(line.split("\t")):
            for a, allele in enumerate(gt.replace("|", "/").split("/")[:ploidy]):
                if allele != ".":
                    G[i, j, a] = int(allele)
    return G


def iter_variant_chunks(vcf_path, chrom, samples, ploidy, chunk_size):
    """
    Yield (pos, G) chunks for SNP records on chrom, streamed from bcftools
    query. pos is 1-based, G is int8 (nsites, nsamples, ploidy) with -1 for
    missing alleles. Records with any non-single-base allele are skipped.
    """
    cmd = ["bcftools", "query", "-r", chrom, "-s", ",".join(samples),
           "-f", "%POS\t%REF,%ALT[\t%GT]\n", vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1 << 20)
    pos_buf, gt_buf = [], []
    for line in proc.stdout:
        pos, alleles, gts = line.rstrip("\n").split("\t", 2)
        n = len(alleles)
        if n % 2 == 0 or alleles[1::2] != "," * (n // 2) or "*" in alleles or "." in alleles:
            continue
        pos_buf.append(int(pos))
        gt_buf.append(gts)
        if len(pos_buf) >= chunk_size:
            yield np.asarray(pos_buf, dtype=np.int64), parse_gt_chunk(gt_buf, len(samples), ploidy)
            pos_buf, gt_buf = [], []
    if pos_buf:
        yield np.asarray(pos_buf, dtype=np.int64), parse_gt_chunk(gt_buf, len(samples), ploidy)
    if proc.wait() != 0:
        sys.exit(f"ERROR: bcftools query failed on {chrom} (is the VCF indexed?)")


def in_mask(pos0, starts, ends):
    """Boolean array: is each 0-based position inside the merged intervals?"""
    if len(starts) == 0:
        return np.zeros(len(pos0), dtype=bool)
    idx = np.searchsorted(starts, pos0, side="right") - 1
    ok = idx >= 0
    res = np.zeros(len(pos0), dtype=bool)
    res[ok] = pos0[ok] < ends[idx[ok]]
    return res


# ----------------------------------------------------------------------------
# Main computation
# ----------------------------------------------------------------------------

def fmt(x):
    return "NA" if x is None or not np.isfinite(x) else f"{x:.6g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True,
                    help="SNP VCF/BCF (bgzipped + indexed)")
    ap.add_argument("--masks", required=True,
                    help="TSV: sample<TAB>mask.bed (callable regions)")
    ap.add_argument("--pops", required=True,
                    help="TSV: sample<TAB>population")
    ap.add_argument("--window-size", type=int, default=100_000,
                    help="window size in bp [100000]")
    chrom_group = ap.add_mutually_exclusive_group()
    chrom_group.add_argument("--chroms", metavar="CHROM[,CHROM...]",
        help=(
            "Comma-separated chromosomes to analyse "
            "[default: all chromosomes represented in masks]"
        ),
    )
    chrom_group.add_argument("--chroms-file", metavar="FILE",
        help=(
            "File with one chromosome name per row "
            "(blank lines and full-line # comments ignored)"
        ),
    )
    ap.add_argument("--ploidy", type=int, default=2,
                    help="ploidy assumed for all samples [2]")
    ap.add_argument("--chunk-size", type=int, default=50_000,
                    help="variant sites per processing chunk")
    ap.add_argument("--stats", nargs="+",
                    choices=["pi", "het", "dxy", "fst", "watterson_theta", "tajima_d"],
                    default=["pi", "het"],
                    help="statistics to estimate [pi het]")
    ap.add_argument("--out", required=True,
                    help="output prefix")
    args = ap.parse_args()

    selected_chroms = None
    try:
        if args.chroms is not None:
            selected_chroms = [name.strip() for name in args.chroms.split(",")]
        elif args.chroms_file is not None:
            selected_chroms = read_chroms_file(args.chroms_file)
    except (OSError, UnicodeError, ValueError) as exc:
        ap.error(f"--chroms-file: {exc}")

    if selected_chroms is not None:
        if not selected_chroms or any(
            not name or any(char.isspace() for char in name)
            for name in selected_chroms
        ):
            ap.error("chromosome selection must contain nonempty names "
                     "without whitespace")
        selected_chroms = list(dict.fromkeys(selected_chroms))

    W, P = args.window_size, args.ploidy
    stats = set(args.stats)
    need_het = bool(stats & {"het", "pi"})
    need_dxy = bool(stats & {"dxy", "fst"})
    need_fst = "fst" in stats
    if need_fst and P != 2:
        sys.exit("ERROR: Weir-Cockerham fst currently requires --ploidy 2")

    # --- samples, populations, masks
    pops_of = dict(read_two_col(args.pops))
    mask_of = dict(read_two_col(args.masks))
    vcf_samples, vcf_contigs = read_vcf_header(args.vcf)

    samples = [s for s in vcf_samples if s in pops_of]
    missing_vcf = [s for s in pops_of if s not in vcf_samples]
    if missing_vcf:
        sys.exit(f"ERROR: samples in pops file not in VCF: {', '.join(missing_vcf)}")
    missing_mask = [s for s in samples if s not in mask_of]
    if missing_mask:
        sys.exit(f"ERROR: samples without a mask: {', '.join(missing_mask)}")

    pop_names = sorted(set(pops_of[s] for s in samples))
    if need_dxy and len(pop_names) < 2:
        sys.exit("ERROR: dxy and fst require at least two populations")
    pop_members = {p: np.array([j for j, s in enumerate(samples) if pops_of[s] == p])
                   for p in pop_names}

    print(f"Loading {len(samples)} masks ...", file=sys.stderr)
    masks = [read_bed(mask_of[s]) for s in samples]

    # --- chromosomes and lengths
    chrom_len = dict(vcf_contigs)
    mask_chroms = []
    for m in masks:
        for c in m:
            if c not in mask_chroms:
                mask_chroms.append(c)
    for c in mask_chroms:
        if chrom_len.get(c) is None:
            chrom_len[c] = int(max(m[c][1][-1] for m in masks if c in m and len(m[c][1])))
    if selected_chroms is not None:
        chroms = selected_chroms
    else:
        order = list(vcf_contigs) + [c for c in mask_chroms if c not in vcf_contigs]
        chroms = [c for c in order if c in set(mask_chroms)]
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64))

    # --- genome-wide accumulators
    gw_diffs = {p: 0.0 for p in pop_names}
    gw_comp = {p: 0.0 for p in pop_names}
    gw_tajima_pi = {p: 0.0 for p in pop_names}
    gw_theta_raw = {p: 0.0 for p in pop_names}
    gw_mutations = {p: 0.0 for p in pop_names}
    gw_callable_sites = {p: 0.0 for p in pop_names}
    gw_tajima_nsum = {p: 0.0 for p in pop_names}
    gw_het = np.zeros(len(samples)) if need_het else None
    gw_call = np.zeros(len(samples)) if need_het else None
    pairs = list(itertools.combinations(pop_names, 2)) if need_dxy else []
    gw_dxy_diffs = {pair: 0.0 for pair in pairs}
    gw_dxy_comp = {pair: 0.0 for pair in pairs}
    gw_wc_a = {pair: 0.0 for pair in pairs}
    gw_wc_b = {pair: 0.0 for pair in pairs}
    gw_wc_c = {pair: 0.0 for pair in pairs}

    fpi = open(f"{args.out}.pi_windows.tsv", "w") if "pi" in stats else None
    fhet = open(f"{args.out}.het_windows.tsv", "w") if "het" in stats else None
    fdxy = open(f"{args.out}.dxy_windows.tsv", "w") if "dxy" in stats else None
    ffst = open(f"{args.out}.fst_windows.tsv", "w") if "fst" in stats else None
    ftheta = open(f"{args.out}.watterson_theta_windows.tsv", "w") if "watterson_theta" in stats else None
    ftajima = open(f"{args.out}.tajima_d_windows.tsv", "w") if "tajima_d" in stats else None
    if fpi:
        fpi.write("chrom\tstart\tend\tpop\tn_snps\tcallable_sites_any\tdiffs\tcomparisons\tpi\n")
    if fhet:
        fhet.write("chrom\tstart\tend\tsample\tpop\tn_het\tn_callable\thet\n")
    if fdxy:
        fdxy.write("chrom\tstart\tend\tpop1\tpop2\tdiffs\tcomparisons\tdxy\n")
    if ffst:
        ffst.write("chrom\tstart\tend\tpop1\tpop2\tfst\n")
    if ftheta:
        ftheta.write("chrom\tstart\tend\tpop\tmutations\tcallable_sites\twatterson_theta\n")
    if ftajima:
        ftajima.write("chrom\tstart\tend\tpop\tmutations\tmean_pairwise_differences\t"
                      "raw_watterson_theta\tobserved_alleles\tcallable_sites\ttajima_d\n")

    for chrom in chroms:
        L = chrom_len.get(chrom)
        if L is None:
            print(f"WARNING: no length for {chrom}, skipping", file=sys.stderr)
            continue
        nwin = (L + W - 1) // W
        print(f"{chrom}: {L} bp, {nwin} windows", file=sys.stderr)

        ivs = [clip(*m.get(chrom, empty), L) for m in masks]

        # callable bases per sample per window
        callable_s = np.vstack([weighted_window_lengths(s, e, np.ones(len(s)), W, nwin)
                                for s, e in ivs])
        het_s = np.zeros_like(callable_s)

        # mask-based haplotype comparisons per population per window
        comp = {}
        any_callable = {}
        for p in pop_names:
            ss, se, sn = coverage_segments([ivs[j] for j in pop_members[p]])
            comp[p] = weighted_window_lengths(ss, se, choose2(P * sn), W, nwin)
            any_callable[p] = weighted_window_lengths(ss, se, np.ones(len(ss)), W, nwin)
        diffs = {p: np.zeros(nwin) for p in pop_names}
        nsnps = {p: np.zeros(nwin) for p in pop_names}
        tajima_pi = {p: np.zeros(nwin) for p in pop_names}
        theta_raw = {p: np.zeros(nwin) for p in pop_names}
        mutations = {p: np.zeros(nwin) for p in pop_names}
        tajima_sites = {p: any_callable[p].copy() for p in pop_names}
        tajima_nsum = {
            p: P * callable_s[pop_members[p]].sum(axis=0)
            for p in pop_names
        }
        dxy_diffs = {pair: np.zeros(nwin) for pair in pairs}
        dxy_comp = {pair: np.zeros(nwin) for pair in pairs}
        wc_a = {pair: np.zeros(nwin) for pair in pairs}
        wc_b = {pair: np.zeros(nwin) for pair in pairs}
        wc_c = {pair: np.zeros(nwin) for pair in pairs}
        for p1, p2 in pairs:
            for j in pop_members[p1]:
                for k in pop_members[p2]:
                    ss, se = intersect_intervals(*ivs[j], *ivs[k])
                    dxy_comp[(p1, p2)] += weighted_window_lengths(
                        ss, se, np.full(len(ss), P * P), W, nwin)

        for pos, G in iter_variant_chunks(args.vcf, chrom, samples, P, args.chunk_size):
            pos0 = pos - 1
            win = pos0 // W
            inm = np.column_stack([in_mask(pos0, s, e) for s, e in ivs])  # (k, nsamp)
            G[~inm] = -1

            if need_het:
                full_call = (G >= 0).all(axis=2)
                het = full_call & (G != G[:, :, :1]).any(axis=2)
                miss_in_mask = inm & ~full_call
                for j in range(len(samples)):
                    het_s[j] += np.bincount(win, weights=het[:, j], minlength=nwin)
                    callable_s[j] -= np.bincount(win, weights=miss_in_mask[:, j], minlength=nwin)

            max_allele = max(int(G.max()), 0) if G.size else 0
            for p in pop_names:
                m = pop_members[p]
                Gp = G[:, m, :].reshape(len(pos), -1)
                h = (Gp >= 0).sum(axis=1)
                same = np.zeros(len(pos))
                allele_counts = []
                for a in range(max_allele + 1):
                    count = (Gp == a).sum(axis=1)
                    allele_counts.append(count)
                    same += choose2(count)
                d = choose2(h) - same
                corr = choose2(h) - choose2(P * inm[:, m].sum(axis=1))
                diffs[p] += np.bincount(win, weights=d, minlength=nwin)
                comp[p] += np.bincount(win, weights=corr, minlength=nwin)
                nsnps[p] += np.bincount(win, weights=d > 0, minlength=nwin)
                site_comp = choose2(h)
                site_pi = np.divide(d, site_comp, out=np.zeros_like(d), where=site_comp > 0)
                tajima_pi[p] += np.bincount(win, weights=site_pi, minlength=nwin)
                allele_count_matrix = np.column_stack(allele_counts)
                n_observed_alleles = (allele_count_matrix > 0).sum(axis=1)
                site_mutations = np.maximum(n_observed_alleles - 1, 0)
                site_theta = np.zeros(len(pos), dtype=np.float64)
                for observed_haps in np.unique(h[site_mutations > 0]):
                    a1, _, _ = tajima_constants(int(observed_haps))
                    if a1:
                        use = (h == observed_haps) & (site_mutations > 0)
                        site_theta[use] = site_mutations[use] / a1
                theta_raw[p] += np.bincount(win, weights=site_theta, minlength=nwin)
                mutations[p] += np.bincount(win, weights=site_mutations, minlength=nwin)
                expected_haps = P * inm[:, m].sum(axis=1)
                tajima_nsum[p] += np.bincount(
                    win, weights=h - expected_haps, minlength=nwin)
                no_called_site = (h == 0) & inm[:, m].any(axis=1)
                tajima_sites[p] -= np.bincount(win, weights=no_called_site, minlength=nwin)

            for p1, p2 in pairs:
                G1 = G[:, pop_members[p1], :].reshape(len(pos), -1)
                G2 = G[:, pop_members[p2], :].reshape(len(pos), -1)
                h1, h2 = (G1 >= 0).sum(axis=1), (G2 >= 0).sum(axis=1)
                cross_comp = h1 * h2
                same = np.zeros(len(pos))
                for a in range(max_allele + 1):
                    same += (G1 == a).sum(axis=1) * (G2 == a).sum(axis=1)
                cross_diffs = cross_comp - same
                mask_comp = P * P * inm[:, pop_members[p1]].sum(axis=1) * inm[:, pop_members[p2]].sum(axis=1)
                dxy_diffs[(p1, p2)] += np.bincount(win, weights=cross_diffs, minlength=nwin)
                dxy_comp[(p1, p2)] += np.bincount(
                    win, weights=cross_comp - mask_comp, minlength=nwin)
                if "fst" in stats:
                    a_site, b_site, c_site, _ = wc_fst_components(
                        G, pop_members[p1], pop_members[p2], max_allele)
                    wc_a[(p1, p2)] += np.bincount(
                        win, weights=a_site, minlength=nwin)
                    wc_b[(p1, p2)] += np.bincount(
                        win, weights=b_site, minlength=nwin)
                    wc_c[(p1, p2)] += np.bincount(
                        win, weights=c_site, minlength=nwin)

        # --- write windows
        for w in range(nwin):
            ws, we = w * W, min((w + 1) * W, L)
            for p in pop_names:
                c = comp[p][w]
                pi = diffs[p][w] / c if c > 0 else None
                if fpi:
                    fpi.write(f"{chrom}\t{ws}\t{we}\t{p}\t{int(nsnps[p][w])}\t"
                              f"{int(any_callable[p][w])}\t{diffs[p][w]:.0f}\t{c:.0f}\t{fmt(pi)}\n")
                if ftheta:
                    theta = theta_raw[p][w] / tajima_sites[p][w] if tajima_sites[p][w] > 0 else None
                    ftheta.write(f"{chrom}\t{ws}\t{we}\t{p}\t{mutations[p][w]:.6g}\t"
                                 f"{int(tajima_sites[p][w])}\t{fmt(theta)}\n")
                if ftajima:
                    stdev = tajima_stdev(tajima_nsum[p][w], tajima_sites[p][w], mutations[p][w])
                    tajima = ((tajima_pi[p][w] - theta_raw[p][w]) / stdev
                              if stdev else None)
                    ftajima.write(f"{chrom}\t{ws}\t{we}\t{p}\t{mutations[p][w]:.6g}\t"
                                  f"{tajima_pi[p][w]:.6g}\t{theta_raw[p][w]:.6g}\t"
                                  f"{tajima_nsum[p][w]:.0f}\t{tajima_sites[p][w]:.0f}\t"
                                  f"{fmt(tajima)}\n")
            for j, s in enumerate(samples):
                if fhet:
                    n = callable_s[j, w]
                    h = het_s[j, w] / n if n > 0 else None
                    fhet.write(f"{chrom}\t{ws}\t{we}\t{s}\t{pops_of[s]}\t"
                               f"{int(het_s[j, w])}\t{int(n)}\t{fmt(h)}\n")
            for p1, p2 in pairs:
                c = dxy_comp[(p1, p2)][w]
                dxy = dxy_diffs[(p1, p2)][w] / c if c > 0 else None
                wc_den = wc_a[(p1, p2)][w] + wc_b[(p1, p2)][w] + wc_c[(p1, p2)][w]
                fst = wc_a[(p1, p2)][w] / wc_den if wc_den > 0 else None
                if fdxy:
                    fdxy.write(f"{chrom}\t{ws}\t{we}\t{p1}\t{p2}\t"
                               f"{dxy_diffs[(p1, p2)][w]:.0f}\t{c:.0f}\t{fmt(dxy)}\n")
                if ffst:
                    ffst.write(f"{chrom}\t{ws}\t{we}\t{p1}\t{p2}\t{fmt(fst)}\n")

        for p in pop_names:
            gw_diffs[p] += diffs[p].sum()
            gw_comp[p] += comp[p].sum()
            gw_tajima_pi[p] += tajima_pi[p].sum()
            gw_theta_raw[p] += theta_raw[p].sum()
            gw_mutations[p] += mutations[p].sum()
            gw_callable_sites[p] += tajima_sites[p].sum()
            gw_tajima_nsum[p] += tajima_nsum[p].sum()
        if need_het:
            gw_het += het_s.sum(axis=1)
            gw_call += callable_s.sum(axis=1)
        for pair in pairs:
            gw_dxy_diffs[pair] += dxy_diffs[pair].sum()
            gw_dxy_comp[pair] += dxy_comp[pair].sum()
            gw_wc_a[pair] += wc_a[pair].sum()
            gw_wc_b[pair] += wc_b[pair].sum()
            gw_wc_c[pair] += wc_c[pair].sum()

    for fh in (fpi, fhet, fdxy, ffst, ftheta, ftajima):
        if fh:
            fh.close()

    # --- genome-wide
    het_gw = None
    if need_het:
        het_gw = np.where(gw_call > 0, gw_het / np.maximum(gw_call, 1), np.nan)
    if "het" in stats:
        with open(f"{args.out}.het_genomewide.tsv", "w") as fh:
            fh.write("sample\tpop\tn_het\tn_callable\thet\n")
            for j, s in enumerate(samples):
                fh.write(f"{s}\t{pops_of[s]}\t{int(gw_het[j])}\t{int(gw_call[j])}\t{fmt(het_gw[j])}\n")
    if "pi" in stats:
        with open(f"{args.out}.pi_genomewide.tsv", "w") as fh:
            fh.write("pop\tn_samples\tdiffs\tcomparisons\tpi\tmean_het\n")
            for p in pop_names:
                c = gw_comp[p]
                pi = gw_diffs[p] / c if c > 0 else None
                members_het = het_gw[pop_members[p]]
                mh = np.nanmean(members_het) if np.isfinite(members_het).any() else None
                fh.write(f"{p}\t{len(pop_members[p])}\t{gw_diffs[p]:.0f}\t{c:.0f}\t{fmt(pi)}\t{fmt(mh)}\n")
    if "watterson_theta" in stats:
        with open(f"{args.out}.watterson_theta_genomewide.tsv", "w") as fh:
            fh.write("pop\tmutations\tcallable_sites\twatterson_theta\n")
            for p in pop_names:
                theta = gw_theta_raw[p] / gw_callable_sites[p] if gw_callable_sites[p] > 0 else None
                fh.write(f"{p}\t{gw_mutations[p]:.6g}\t{int(gw_callable_sites[p])}\t{fmt(theta)}\n")
    if "tajima_d" in stats:
        with open(f"{args.out}.tajima_d_genomewide.tsv", "w") as fh:
            fh.write("pop\tmutations\tmean_pairwise_differences\traw_watterson_theta\t"
                     "observed_alleles\tcallable_sites\ttajima_d\n")
            for p in pop_names:
                stdev = tajima_stdev(gw_tajima_nsum[p], gw_callable_sites[p], gw_mutations[p])
                tajima = ((gw_tajima_pi[p] - gw_theta_raw[p]) / stdev
                          if stdev else None)
                fh.write(f"{p}\t{gw_mutations[p]:.6g}\t{gw_tajima_pi[p]:.6g}\t"
                         f"{gw_theta_raw[p]:.6g}\t{gw_tajima_nsum[p]:.0f}\t"
                         f"{gw_callable_sites[p]:.0f}\t{fmt(tajima)}\n")
    if "dxy" in stats:
        with open(f"{args.out}.dxy_genomewide.tsv", "w") as fh:
            fh.write("pop1\tpop2\tdiffs\tcomparisons\tdxy\n")
            for p1, p2 in pairs:
                c = gw_dxy_comp[(p1, p2)]
                dxy = gw_dxy_diffs[(p1, p2)] / c if c > 0 else None
                fh.write(f"{p1}\t{p2}\t{gw_dxy_diffs[(p1, p2)]:.0f}\t{c:.0f}\t{fmt(dxy)}\n")
    if "fst" in stats:
        with open(f"{args.out}.fst_genomewide.tsv", "w") as fh:
            fh.write("pop1\tpop2\tfst\n")
            for p1, p2 in pairs:
                c = gw_dxy_comp[(p1, p2)]
                dxy = gw_dxy_diffs[(p1, p2)] / c if c > 0 else None
                wc_den = gw_wc_a[(p1, p2)] + gw_wc_b[(p1, p2)] + gw_wc_c[(p1, p2)]
                fst = gw_wc_a[(p1, p2)] / wc_den if wc_den > 0 else None
                fh.write(f"{p1}\t{p2}\t{fmt(fst)}\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
