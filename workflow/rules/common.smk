import math
import os
import re

import pandas as pd

SCRIPTS = os.path.join(workflow.basedir, "scripts")
OUT = config.get("outdir", "results")
ENV_DIR = f"{OUT}/env"

SAMPLES = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .fillna("")
    .set_index("species", drop=False)
)
for col in ("species", "vcf", "popmap", "latlongs"):
    if col not in SAMPLES.columns:
        raise ValueError(f"sample sheet {config['samples']} is missing column '{col}'")
if "id_regex" not in SAMPLES.columns:
    SAMPLES["id_regex"] = ""
bad = [s for s in SAMPLES.index if not re.fullmatch(r"[A-Za-z0-9_.-]+", s)]
if bad:
    raise ValueError(f"species names may only contain letters, digits, '_', '.', '-': {bad}")
if SAMPLES.index.duplicated().any():
    raise ValueError("duplicate species names in sample sheet")

SPECIES = list(SAMPLES.index)
N_RANDOM = int(config["gf"]["n_random"])
RANDOM_BATCH = max(1, int(config["gf"].get("random_batch_size", 10)))
N_BATCHES = -(-N_RANDOM // RANDOM_BATCH)  # ceiling


def random_runtime(wc, input):
    """SLURM minutes for one randomization batch: ~3 min per 100 trees per fit
    (about 14k SNPs, 14 sites) doubled for margin, times the batch size."""
    ntree = int(pd.read_csv(input.ntree, sep="\t")["ntree"].iloc[0])
    return 30 + RANDOM_BATCH * math.ceil(6 * ntree / 100)


def batch_reps(wc):
    """first and last replicate number in a randomization batch"""
    b = int(wc.batch)
    return (b - 1) * RANDOM_BATCH + 1, min(b * RANDOM_BATCH, N_RANDOM)


wildcard_constraints:
    species="|".join(re.escape(s) for s in SPECIES),
    batch=r"\d+",


def sample_col(col):
    return lambda wc: SAMPLES.loc[wc.species, col]


def sp_dir(sub=""):
    return f"{OUT}/{{species}}/{sub}".rstrip("/")


def all_outputs():
    plots = [
        "map_turnover_rgb",
        "biplot_turnover_pca",
        "importance_r2weighted",
        "importance_accuracy",
        "random_vs_real",
    ]
    out = []
    for sp in SPECIES:
        for p in plots:
            out += [f"{OUT}/{sp}/plots/{p}.pdf", f"{OUT}/{sp}/plots/{p}.png"]
        out += [
            f"{OUT}/{sp}/plots/gf_diagnostics.pdf",
            f"{OUT}/{sp}/qc.html",
            f"{OUT}/{sp}/gf/random_pvalues.tsv",
            f"{OUT}/{sp}/predictors/correlation.png",
        ]
    return out
