#!/usr/bin/env python
"""PCA of individual genotypes (population structure and QC; not used by gradient forest).

Individuals: every VCF sample that maps to a population with coordinates
(qc/sample_pop.tsv), including populations later dropped from the gradient
forest model for having too few individuals (column in_gf_model).
SNPs: biallelic SNPs whose genotypes vary among those individuals. No MAF or missingness
filter.

Each SNP is centred on twice its allele frequency and scaled by
sqrt(2p(1 - p)) (Patterson et al. 2006, as in EIGENSOFT and PLINK); missing
genotypes are set to the mean (0 after centring). The PCs are the eigenvectors
of the individuals x individuals matrix X X' / n_snps, accumulated over the VCF
in blocks of SNPs, so memory depends on the number of individuals, not SNPs.

Outputs
  --out-variance  PC, eigenvalue, pct_variance, cumulative_pct
  --out-scores    one row per individual: sample_id, popid, sitename,
                  missing_frac, in_gf_model, PC1..PCk (coordinates on each PC,
                  = eigenvector x sqrt(eigenvalue)); signs are fixed so the
                  individual farthest from 0 on each PC is positive
  --out-info      counts of individuals and SNPs
Individuals with much missing data are pulled towards 0 on every PC; check
missing_frac before reading structure into such samples.
"""
import argparse
import sys
from collections import Counter

import numpy as np
import pandas as pd
from cyvcf2 import VCF

BLOCK = 5000  # SNPs per block added to X X'


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--sample-pop", required=True)
    p.add_argument("--env-sites", required=True)
    p.add_argument("--pop-counts", required=True, help="genotypes/pop_counts.tsv (populations in the GF model)")
    p.add_argument("--out-variance", required=True)
    p.add_argument("--out-scores", required=True)
    p.add_argument("--out-info", required=True)
    a = p.parse_args()

    sp = pd.read_csv(a.sample_pop, sep="\t", dtype=str)
    sites = pd.read_csv(a.env_sites, sep="\t", dtype={"popid": str, "sitename": str})
    gf_pops = set(pd.read_csv(a.pop_counts, sep="\t", dtype={"popid": str}).popid)
    if len(sp) < 3:
        sys.exit("ERROR: fewer than 3 individuals; no PCA")

    vcf = VCF(a.vcf, samples=list(sp.vcf_id), gts012=False, lazy=True, threads=2)
    samples = list(vcf.samples)  # VCF order
    n = len(samples)
    K = np.zeros((n, n))
    n_missing = np.zeros(n)
    counts = Counter()
    block = []

    def add_block():
        nonlocal K
        X = np.column_stack(block)
        K += X @ X.T
        block.clear()

    for v in vcf:
        counts["total"] += 1
        if len(v.ALT) != 1 or not v.is_snp:
            counts["not_biallelic_snp"] += 1
            continue
        gt = v.gt_types  # 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
        called = gt != 2
        if not called.any():
            counts["no_calls"] += 1
            continue
        g = np.where(gt == 3, 2.0, np.where(gt == 1, 1.0, 0.0))
        if np.ptp(g[called]) == 0:  # every called individual has the same genotype
            counts["invariant"] += 1
            continue
        f = g[called].mean() / 2
        x = np.where(called, (g - 2 * f) / np.sqrt(2 * f * (1 - f)), 0.0)
        block.append(x)
        n_missing += ~called
        counts["used"] += 1
        if len(block) == BLOCK:
            add_block()
    if block:
        add_block()
    m = counts["used"]
    if m < 2:
        sys.exit(f"ERROR: only {m} variable SNPs; no PCA")

    lam, U = np.linalg.eigh(K / m)
    order = np.argsort(lam)[::-1]
    lam, U = lam[order], U[:, order]
    keep = lam > lam[0] * 1e-10  # centring leaves at most n - 1 non-zero eigenvalues
    lam, U = lam[keep], U[:, keep]
    scores = U * np.sqrt(lam)
    flip = np.sign(scores[np.abs(scores).argmax(axis=0), np.arange(scores.shape[1])])
    scores *= flip
    k = len(lam)
    pcs = [f"PC{i + 1}" for i in range(k)]

    total = np.trace(K / m)  # total variance = sum of all eigenvalues
    var = pd.DataFrame({"PC": pcs, "eigenvalue": lam, "pct_variance": 100 * lam / total})
    var["cumulative_pct"] = var.pct_variance.cumsum()
    var.to_csv(a.out_variance, sep="\t", index=False, float_format="%.6g")

    pop_of = dict(zip(sp.vcf_id, sp.popid))
    site_of = dict(zip(sites.popid, sites.sitename))
    pops = [pop_of[s] for s in samples]
    out = pd.DataFrame({
        "sample_id": samples,
        "popid": pops,
        "sitename": [site_of.get(x, "") for x in pops],
        "missing_frac": n_missing / m,
        "in_gf_model": ["yes" if x in gf_pops else "no" for x in pops],
    })
    out = pd.concat([out, pd.DataFrame(scores, columns=pcs)], axis=1)
    out.to_csv(a.out_scores, sep="\t", index=False, float_format="%.6g")

    names = ["total", "not_biallelic_snp", "no_calls", "invariant", "used"]
    info = pd.DataFrame({"item": ["n_individuals"] + [f"snps_{c}" for c in names],
                         "n": [n] + [counts[c] for c in names]})
    info.to_csv(a.out_info, sep="\t", index=False)
    print(info.to_string(index=False))
    print(var.head(10).to_string(index=False))
    hi = out[out.missing_frac > 0.5]
    if len(hi):
        print(f"WARNING: {len(hi)} individuals with > 50% missing genotypes: "
              f"{', '.join(hi.sample_id[:10])}{' ...' if len(hi) > 10 else ''}", file=sys.stderr)


if __name__ == "__main__":
    main()
