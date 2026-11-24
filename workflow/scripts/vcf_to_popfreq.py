#!/usr/bin/env python
"""Population allele frequencies from a VCF, as gradient forest response variables.

Keeps biallelic SNPs; drops populations with fewer than --min-ind-per-pop
individuals. No MAF or missingness filter: gradient forest itself leaves out
SNPs without enough signal (R^2 <= 0). The only SNPs dropped are those with no
called genotype in some population (frequency undefined) and those with no
frequency variation among populations (nothing to model).
Frequencies are of the globally minor allele.
Output: rows = populations (in env_sites order), columns = SNPs; --out-snp-info
maps each SNP column name to its VCF CHROM, POS and ID.
"""
import argparse
import sys
from collections import Counter

import numpy as np
import pandas as pd
from cyvcf2 import VCF


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--sample-pop", required=True)
    p.add_argument("--env-sites", required=True)
    p.add_argument("--min-ind-per-pop", type=int, default=5)
    p.add_argument("--out", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-pops", required=True)
    p.add_argument("--out-snp-info", required=True)
    a = p.parse_args()

    sp = pd.read_csv(a.sample_pop, sep="\t", dtype=str)
    sites = pd.read_csv(a.env_sites, sep="\t", dtype={"popid": str})
    sizes = sp.popid.value_counts()
    small = [x for x in sites.popid if sizes.get(x, 0) < a.min_ind_per_pop]
    if small:
        print(f"WARNING: dropping populations with < {a.min_ind_per_pop} individuals: {small}", file=sys.stderr)
    pops = [x for x in sites.popid if x not in small]
    sp = sp[sp.popid.isin(pops)]
    if len(pops) < 3:
        sys.exit("ERROR: fewer than 3 populations remain")

    vcf = VCF(a.vcf, samples=list(sp.vcf_id), gts012=False, lazy=True, threads=2)
    pop_of = dict(zip(sp.vcf_id, sp.popid))
    pidx = np.array([pops.index(pop_of[s]) for s in vcf.samples])
    P = len(pops)
    n_ind = np.bincount(pidx, minlength=P).astype(float)

    reasons = Counter()
    ids, freqs, info = [], [], []
    seen = Counter()
    for v in vcf:
        reasons["total"] += 1
        if len(v.ALT) != 1 or not v.is_snp:
            reasons["not_biallelic_snp"] += 1
            continue
        gt = v.gt_types  # 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
        called = gt != 2
        alt = np.where(gt == 1, 1.0, np.where(gt == 3, 2.0, 0.0))
        n_called = np.bincount(pidx, weights=called.astype(float), minlength=P)
        if np.any(n_called == 0):
            reasons["no_calls_in_a_pop"] += 1
            continue
        n_alt = np.bincount(pidx, weights=alt, minlength=P)
        p_glob = n_alt.sum() / (2 * n_called.sum())
        f = n_alt / (2 * n_called)
        if p_glob > 0.5:
            f = 1 - f
        if np.ptp(f) == 0:
            reasons["no_among_pop_variation"] += 1
            continue
        sid = f"snp_{v.ID}" if v.ID not in (None, ".") else f"{v.CHROM}_{v.POS}"
        seen[sid] += 1
        if seen[sid] > 1:
            sid = f"{sid}_{seen[sid]}"
        ids.append(sid)
        info.append((sid, v.CHROM, v.POS, v.ID if v.ID is not None else "."))
        freqs.append(f.astype("float32"))
        reasons["kept"] += 1

    if not freqs:
        sys.exit("ERROR: no SNPs left")
    mat = np.column_stack(freqs)
    df = pd.DataFrame(mat, columns=ids)
    df.insert(0, "popid", pops)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4f")

    order = ["total", "not_biallelic_snp", "no_calls_in_a_pop", "no_among_pop_variation", "kept"]
    summ = pd.DataFrame({"category": order, "n_snps": [reasons[k] for k in order]})
    summ.to_csv(a.out_summary, sep="\t", index=False)
    pd.DataFrame({"popid": pops, "n_ind": n_ind.astype(int)}).to_csv(a.out_pops, sep="\t", index=False)
    print(summ.to_string(index=False))
    pd.DataFrame(info, columns=["snp", "CHROM", "POS", "VCF_ID"]).to_csv(a.out_snp_info, sep="\t", index=False)
    print(f"matrix: {P} populations x {len(ids)} SNPs")


if __name__ == "__main__":
    main()
