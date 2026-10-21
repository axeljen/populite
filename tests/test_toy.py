#!/usr/bin/env python3
"""
Random toy data test: compare populite.py with a brute-force
per-base implementation. Run: python tests/test_toy.py
"""

import os
import random
import subprocess
import sys
import tempfile

import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "populite.py")
SPEC = importlib.util.spec_from_file_location("populite", TOOL)
SUMSTATS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMSTATS)

CHROMS = {"chr1": 2750, "chr2": 1200}
WSIZE = 1000
SAMPLES = [f"s{i}" for i in range(7)]
POPS = {s: ("A" if i < 4 else "B") for i, s in enumerate(SAMPLES)}


def c2(x):
    return x * (x - 1) / 2


def make_data(d, rng):
    # masks: random intervals, some overlapping, some spanning window edges
    masks = {}
    for s in SAMPLES:
        masks[s] = {}
        with open(os.path.join(d, f"{s}.bed"), "w") as fh:
            for c, L in CHROMS.items():
                cov = set()
                for _ in range(rng.randint(3, 8)):
                    a = rng.randrange(0, L)
                    b = min(L + 50, a + rng.randint(1, 900))  # may overhang the end
                    fh.write(f"{c}\t{a}\t{b}\n")
                    cov.update(range(a, min(b, L)))
                masks[s][c] = cov
    # variants: 1-based positions, genotypes as lists of allele ints or None
    variants = {}
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={L}>" for c, L in CHROMS.items()]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES))
    for c, L in CHROMS.items():
        variants[c] = {}
        for pos in sorted(rng.sample(range(1, L + 1), 120)):
            nalt = 2 if rng.random() < 0.1 else 1
            gts, gtstr = [], []
            for _ in SAMPLES:
                if rng.random() < 0.1:
                    gts.append(None)
                    gtstr.append("./.")
                else:
                    g = [rng.randint(0, nalt) for _ in range(2)]
                    gts.append(g)
                    gtstr.append(f"{g[0]}/{g[1]}")
            variants[c][pos] = gts
            alt = "C,G" if nalt == 2 else "C"
            lines.append(f"{c}\t{pos}\t.\tA\t{alt}\t50\tPASS\t.\tGT\t" + "\t".join(gtstr))
        # an indel that must be ignored
    lines.append(f"chr1\t{CHROMS['chr1']}\t.\tAT\tA\t50\tPASS\t.\tGT\t" + "\t".join(["0/1"] * len(SAMPLES)))
    lines = lines[:5] + sorted(lines[5:], key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
    vcf = os.path.join(d, "toy.vcf")
    with open(vcf, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    subprocess.run(["bcftools", "view", "-Oz", "-o", vcf + ".gz", vcf], check=True)
    subprocess.run(["bcftools", "index", "-t", vcf + ".gz"], check=True)
    with open(os.path.join(d, "masks.tsv"), "w") as fh:
        for s in SAMPLES:
            fh.write(f"{s}\t{os.path.join(d, s + '.bed')}\n")
    with open(os.path.join(d, "pops.tsv"), "w") as fh:
        for s in SAMPLES:
            fh.write(f"{s}\t{POPS[s]}\n")
    # the indel at chr1:L is not a SNP; drop it from the truth if present
    return masks, variants


def brute_force(masks, variants):
    pi_win, het_win, dxy_win, theta_win, tajima_win = {}, {}, {}, {}, {}
    for c, L in CHROMS.items():
        for w in range((L + WSIZE - 1) // WSIZE):
            dxy_win[(c, w)] = [0, 0]
            for p in "AB":
                members = [s for s in SAMPLES if POPS[s] == p]
                diffs = comp = 0
                for b in range(w * WSIZE, min((w + 1) * WSIZE, L)):
                    haps = []
                    for s in members:
                        if b not in masks[s][c]:
                            continue
                        gts = variants[c].get(b + 1)
                        if gts is None:
                            haps += [0, 0]
                        else:
                            g = gts[SAMPLES.index(s)]
                            if g is not None:
                                haps += g
                    n = len(haps)
                    comp += c2(n)
                    diffs += c2(n) - sum(c2(haps.count(a)) for a in set(haps))
                pi_win[(c, w, p)] = (diffs, comp)
            for i, s in enumerate(SAMPLES):
                nh = nc = 0
                for b in range(w * WSIZE, min((w + 1) * WSIZE, L)):
                    if b not in masks[s][c]:
                        continue
                    gts = variants[c].get(b + 1)
                    if gts is None:
                        nc += 1
                    elif gts[i] is not None:
                        nc += 1
                        nh += gts[i][0] != gts[i][1]
                het_win[(c, w, s)] = (nh, nc)
            for b in range(w * WSIZE, min((w + 1) * WSIZE, L)):
                gts = variants[c].get(b + 1)
                pop_haps = {}
                for p in "AB":
                    pop_haps[p] = []
                    for i, s in enumerate(SAMPLES):
                        if POPS[s] != p or b not in masks[s][c]:
                            continue
                        g = None if gts is None else gts[i]
                        if gts is None:
                            pop_haps[p] += [0, 0]
                        elif g is not None:
                            pop_haps[p] += g
                h1, h2 = pop_haps["A"], pop_haps["B"]
                dxy_win[(c, w)][0] += sum(a != b for a in h1 for b in h2)
                dxy_win[(c, w)][1] += len(h1) * len(h2)
            for p in "AB":
                members = [i for i, s in enumerate(SAMPLES) if POPS[s] == p]
                raw_pi = raw_theta = total_alleles = sites = mutations = 0
                for b in range(w * WSIZE, min((w + 1) * WSIZE, L)):
                    gts = variants[c].get(b + 1)
                    haps = []
                    for i in members:
                        s = SAMPLES[i]
                        if b not in masks[s][c]:
                            continue
                        if gts is None:
                            haps += [0, 0]
                        elif gts[i] is not None:
                            haps += gts[i]
                    n = len(haps)
                    if n == 0:
                        continue
                    sites += 1
                    total_alleles += n
                    if n > 1:
                        diffs = c2(n) - sum(c2(haps.count(a)) for a in set(haps))
                        raw_pi += diffs / c2(n)
                    n_mutations = max(len(set(haps)) - 1, 0)
                    mutations += n_mutations
                    if n_mutations:
                        a1 = sum(1 / i for i in range(1, n))
                        raw_theta += n_mutations / a1
                n_mean = round(total_alleles / sites) if sites else 0
                _, e1, e2 = SUMSTATS.tajima_constants(n_mean)
                variance = (e1 * mutations + e2 * mutations * (mutations - 1)
                            if e1 is not None and mutations else 0)
                tajima_win[(c, w, p)] = ((raw_pi - raw_theta) / variance ** 0.5
                                         if variance > 0 else None)
                theta_win[(c, w, p)] = raw_theta / sites if sites else None
    return pi_win, het_win, dxy_win, theta_win, tajima_win


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    g = np.array([[[0, 0], [0, 0], [1, 1], [1, 1]]], dtype=np.int8)
    a, b, c, valid = SUMSTATS.wc_fst_components(g, np.array([0, 1]),
                                               np.array([2, 3]), 1)
    if not valid[0] or abs(a[0] / (a[0] + b[0] + c[0]) - 1) > 1e-12:
        print("WC FST FIXED-DIFFERENCE MISMATCH")
        sys.exit(1)
    with tempfile.TemporaryDirectory() as d:
        masks, variants = make_data(d, rng)
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, TOOL, "--vcf", os.path.join(d, "toy.vcf.gz"),
                        "--masks", os.path.join(d, "masks.tsv"),
                        "--pops", os.path.join(d, "pops.tsv"),
                        "--window-size", str(WSIZE), "--chunk-size", "37",
                        "--out", out], check=True)
        pi_true, het_true, dxy_true, theta_true, tajima_true = brute_force(masks, variants)

        nfail = 0
        with open(out + ".pi_windows.tsv") as fh:
            next(fh)
            for line in fh:
                c, s, e, p, _, _, d_, cmp_, _ = line.split("\t")
                exp = pi_true[(c, int(s) // WSIZE, p)]
                if (float(d_), float(cmp_)) != exp:
                    nfail += 1
                    print("PI MISMATCH", c, s, p, d_, cmp_, exp)
        with open(out + ".het_windows.tsv") as fh:
            next(fh)
            for line in fh:
                c, s, e, smp, _, nh, nc, _ = line.split("\t")
                exp = het_true[(c, int(s) // WSIZE, smp)]
                if (int(nh), int(nc)) != exp:
                    nfail += 1
                    print("HET MISMATCH", c, s, smp, nh, nc, exp)
        # genome-wide must equal sums of windows
        gd = {p: sum(v[0] for k, v in pi_true.items() if k[2] == p) for p in "AB"}
        gc = {p: sum(v[1] for k, v in pi_true.items() if k[2] == p) for p in "AB"}
        with open(out + ".pi_genomewide.tsv") as fh:
            next(fh)
            for line in fh:
                p, _, d_, cmp_ = line.split("\t")[:4]
                if (float(d_), float(cmp_)) != (gd[p], gc[p]):
                    nfail += 1
                    print("GW PI MISMATCH", p, d_, cmp_, gd[p], gc[p])

        stats_out = os.path.join(d, "allstats")
        subprocess.run([sys.executable, TOOL, "--vcf", os.path.join(d, "toy.vcf.gz"),
                        "--masks", os.path.join(d, "masks.tsv"),
                        "--pops", os.path.join(d, "pops.tsv"), "--window-size", str(WSIZE),
                        "--chunk-size", "37", "--stats", "pi", "het", "dxy", "fst",
                        "watterson_theta", "tajima_d", "--out", stats_out], check=True)
        with open(stats_out + ".dxy_windows.tsv") as fh:
            next(fh)
            for line in fh:
                c, start, end, p1, p2, diffs, comp, dxy = line.split("\t")
                exp = dxy_true[(c, int(start) // WSIZE)]
                if (int(diffs), int(comp)) != tuple(exp):
                    nfail += 1
                    print("DXY MISMATCH", c, start, diffs, comp, exp)
                expected_dxy = exp[0] / exp[1] if exp[1] else None
                if expected_dxy is None and dxy != "NA":
                    nfail += 1
                elif expected_dxy is not None and abs(float(dxy) - expected_dxy) > 1e-6:
                    nfail += 1
                    print("DXY RATE MISMATCH", c, start, dxy, expected_dxy)
        with open(stats_out + ".watterson_theta_windows.tsv") as fh:
            next(fh)
            for line in fh:
                c, start, end, p, mutations, sites, theta = line.split("\t")
                expected = theta_true[(c, int(start) // WSIZE, p)]
                if expected is not None and abs(float(theta) - expected) > 1e-6:
                    nfail += 1
                    print("THETA MISMATCH", c, start, p, theta, expected)
        with open(stats_out + ".tajima_d_windows.tsv") as fh:
            next(fh)
            for line in fh:
                c, start, end, p, mutations, raw_pi, raw_theta, nsum, sites, tajima = line.strip().split("\t")
                expected = tajima_true[(c, int(start) // WSIZE, p)]
                if expected is None and tajima != "NA":
                    nfail += 1
                    print("TAJIMA D UNEXPECTED VALUE", c, start, p, tajima)
                elif expected is not None and abs(float(tajima) - expected) > 1e-5:
                    nfail += 1
                    print("TAJIMA D MISMATCH", c, start, p, tajima, expected)
        theta_only = os.path.join(d, "theta_only")
        subprocess.run([sys.executable, TOOL, "--vcf", os.path.join(d, "toy.vcf.gz"),
                        "--masks", os.path.join(d, "masks.tsv"),
                        "--pops", os.path.join(d, "pops.tsv"), "--stats", "watterson_theta",
                        "--out", theta_only], check=True)
        if (os.path.exists(theta_only + ".pi_windows.tsv") or
                os.path.exists(theta_only + ".het_windows.tsv") or
                not os.path.exists(theta_only + ".watterson_theta_genomewide.tsv")):
            nfail += 1
            print("STAT SELECTION MISMATCH")

        single_pops = os.path.join(d, "single_pops.tsv")
        with open(single_pops, "w") as fh:
            for s in SAMPLES:
                fh.write(f"{s}\tA\n")
        result = subprocess.run([sys.executable, TOOL, "--vcf", os.path.join(d, "toy.vcf.gz"),
                                 "--masks", os.path.join(d, "masks.tsv"),
                                 "--pops", single_pops, "--stats", "fst", "--out",
                                 os.path.join(d, "single")], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode == 0 or "require at least two populations" not in result.stderr:
            nfail += 1
            print("SINGLE-POPULATION FST VALIDATION MISMATCH")
        print("FAIL" if nfail else "PASS", f"({nfail} mismatches)")
        sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
