#!/usr/bin/env python
"""Generate the small synthetic dataset in .test/data used by the CI tests.

The generated files are committed, so this only needs rerunning when the test
data have to change. Needs only numpy. Run from the repository root:

    python .test/make_test_data.py

Study area: a 4 x 3 degree box (lon 10-14, lat 0-3) on a 0.1-degree grid.
Six environmental layers, written as ESRI ASCII grids:
  bio1     east-west gradient (drives most of the signal in species A)
  bio_2    bio1 plus a little noise, so |r| with bio1 > 0.8 at the sites
           (file name with an underscore, to test layer-name normalization)
  bio3     north-south gradient (drives species B)
  bio4     smooth wave pattern (drives part of species B)
  bio5     constant everywhere (must be dropped as "constant across sites")
  bio6     spatial noise (no signal)
The north-west corner (lon < 10.5, lat > 2.5) is NODATA in every layer; one
species A site sits just inside it and must be snapped to the nearest valid cell.

Species A (spA): 13 populations in the latlongs file, 10 individuals each,
except population 13 with 3 (dropped by min_ind_per_pop = 4). VCF sample IDs
carry a "TEST_" prefix removed by id_regex; one extra VCF sample is not in the
popmap (must be listed as unmapped). Plain-text VCF, 400 SNPs plus one
multiallelic site and one SNP with no variation among populations.

Species B (spB): 10 populations of 8 individuals, gzip-compressed VCF, 250 SNPs.
"""
import gzip
import os

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
XLL, YLL, CELL, NCOLS, NROWS = 10.0, 0.0, 0.1, 40, 30
NODATA = -9999


def grid_xy():
    xs = XLL + CELL * (np.arange(NCOLS) + 0.5)
    ys = YLL + CELL * (NROWS - np.arange(NROWS) - 0.5)  # first row = north
    return np.meshgrid(xs, ys)


def layers(rng):
    x, y = grid_xy()
    bio1 = 20 + 2.0 * (x - XLL)
    return {
        "bio1": bio1,
        "bio_2": 0.9 * bio1 + rng.normal(0, 0.3, x.shape),
        "bio3": 100 + 40 * (y - YLL),
        "bio4": 50 + 10 * np.sin(1.7 * x) * np.cos(2.3 * y),
        "bio5": np.full(x.shape, 7.0),
        "bio6": rng.normal(0, 1, x.shape),
    }


def write_asc(path, arr, mask):
    a = np.where(mask, NODATA, arr)
    with open(path, "w") as fh:
        fh.write(f"ncols {NCOLS}\nnrows {NROWS}\nxllcorner {XLL}\nyllcorner {YLL}\n"
                 f"cellsize {CELL}\nNODATA_value {NODATA}\n")
        for row in a:
            fh.write(" ".join(f"{v:.4f}" if v != NODATA else str(NODATA) for v in row) + "\n")


def env_at(lyr, lon, lat):
    c = int((lon - XLL) / CELL)
    r = int((YLL + NROWS * CELL - lat) / CELL)
    return {k: v[r, c] for k, v in lyr.items()}


def logistic(z):
    return 1 / (1 + np.exp(-z))


def scale(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std()


def pop_freqs(rng, env, n_snp, drivers):
    """allele frequencies (pops x SNPs): a share of SNPs follow a cline along
    each driver variable; the rest are random among populations"""
    n_pop = len(next(iter(env.values())))
    f = rng.uniform(0.1, 0.9, (n_pop, n_snp))
    start = 0
    for var, share in drivers:
        k = int(n_snp * share)
        z = scale(env[var])
        slope = rng.choice([-1, 1], k) * rng.uniform(1.5, 3.0, k)
        shift = rng.normal(0, 0.3, k)
        f[:, start:start + k] = 0.05 + 0.9 * logistic(np.outer(z, slope) + shift)
        start += k
    return f


def genotype_lines(rng, freqs, pops_of_samples, chrom, missing=0.05):
    """VCF body lines: one SNP per column of freqs"""
    lines = []
    for j in range(freqs.shape[1]):
        p = freqs[pops_of_samples, j]
        g = rng.binomial(2, p)
        gts = np.array(["0/0", "0/1", "1/1"])[g]
        gts[rng.random(len(g)) < missing] = "./."
        lines.append(f"{chrom}\t{1000 * (j + 1)}\t{chrom}_{j + 1}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(gts))
    return lines


def vcf_header(samples):
    return ["##fileformat=VCFv4.2", "##contig=<ID=chr1>", "##contig=<ID=chr2>",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]


def make_species(rng, lyr, name, sites, n_ind, n_snp, drivers, prefix, chrom,
                 extra_samples=(), extra_lines=()):
    """latlongs and popmap files; returns (directory, VCF sample names, VCF body lines).
    extra_samples are in the VCF but not the popmap; extra_lines are VCF sites
    (columns up to FORMAT) given 0/1 in every sample."""
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{name}_pop_latlongs.csv"), "w") as fh:
        fh.write("popid,sitename,lat,long\n")
        for i, (lon, lat) in enumerate(sites, start=1):
            fh.write(f"{i},{name}_site{i},{lat:.4f},{lon:.4f}\n")
    # env at each site: the snapped site uses its nearest valid cell (one column east)
    env = {k: [] for k in lyr}
    for lon, lat in sites:
        e = env_at(lyr, lon, lat)
        if np.isnan(e["bio1"]):
            e = env_at(lyr, lon + CELL, lat)
        for k in lyr:
            env[k].append(e[k])
    freqs = pop_freqs(rng, env, n_snp, drivers)

    samples, pops_idx, popmap = [], [], []
    for i, n in enumerate(n_ind):
        for s in range(1, n + 1):
            sid = f"{name}_p{i + 1:02d}_{s:02d}"
            samples.append(prefix + sid)
            pops_idx.append(i)
            popmap.append((sid, i + 1))
    body = genotype_lines(rng, freqs, np.array(pops_idx), chrom)
    samples += list(extra_samples)
    body = [line + "\t0/1" * len(extra_samples) for line in body]
    body += [line + "\t0/1" * len(samples) for line in extra_lines]
    with open(os.path.join(d, f"{name}_id_to_pop.txt"), "w") as fh:
        fh.write("sampleid popid\n")
        fh.writelines(f"{s} {p}\n" for s, p in popmap)
    return d, samples, body


def main():
    rng = np.random.default_rng(20261125)
    lyr = layers(rng)
    x, y = grid_xy()
    mask = (x < 10.5) & (y > 2.5)
    for k in lyr:
        lyr[k] = np.where(mask, np.nan, lyr[k])
    env_dir = os.path.join(OUT, "env")
    os.makedirs(env_dir, exist_ok=True)
    for k, v in lyr.items():
        write_asc(os.path.join(env_dir, f"test-{k}.asc"), v, mask)
    with open(os.path.join(OUT, "study_area.geojson"), "w") as fh:
        fh.write('{"type": "FeatureCollection", "features": [{"type": "Feature", "properties": '
                 '{"name": "test"}, "geometry": {"type": "Polygon", "coordinates": '
                 '[[[10.0, 0.0], [14.0, 0.0], [14.0, 3.0], [10.0, 3.0], [10.0, 0.0]]]}}]}\n')

    # species A: 12 sites on a jittered grid, the 13th (3 individuals) elsewhere;
    # site 1 lies in the NODATA corner, 0.02 degrees from valid cells
    sites_a = [(10.48, 2.62)]
    for lon in (11.0, 12.0, 13.0, 13.5):
        for lat in (0.5, 1.5, 2.5):
            if len(sites_a) < 12:
                sites_a.append((lon + rng.uniform(-0.3, 0.3), lat + rng.uniform(-0.3, 0.3)))
    sites_a.append((12.5, 2.0))
    n_ind_a = [10] * 12 + [3]
    extra = ["chr1\t999001\tmulti_1\tA\tG,T\t.\tPASS\t.\tGT",   # multiallelic: dropped
             "chr1\t999002\tflat_1\tC\tT\t.\tPASS\t.\tGT"]      # no variation among pops: dropped
    d, samples, body = make_species(rng, lyr, "spA", sites_a, n_ind_a, 400,
                                    [("bio1", 0.35), ("bio3", 0.15)], "TEST_", "chr1",
                                    extra_samples=["TEST_spA_extra_01"], extra_lines=extra)
    with open(os.path.join(d, "spA.vcf"), "w") as fh:
        fh.write("\n".join(vcf_header(samples) + body) + "\n")

    # species B: 10 sites, 8 individuals each
    sites_b = [(rng.uniform(10.8, 13.8), rng.uniform(0.3, 2.7)) for _ in range(10)]
    d, samples, body = make_species(rng, lyr, "spB", sites_b, [8] * 10, 250,
                                    [("bio3", 0.3), ("bio4", 0.2)], "", "chr2")
    with gzip.open(os.path.join(d, "spB.vcf.gz"), "wt") as fh:
        fh.write("\n".join(vcf_header(samples) + body) + "\n")
    print(f"wrote test data to {OUT}")


if __name__ == "__main__":
    main()
