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
for col in ("id_regex", "min_ind_per_pop"):  # optional columns
    if col not in SAMPLES.columns:
        SAMPLES[col] = ""
bad_n = [s for s, v in SAMPLES.min_ind_per_pop.items() if v and not (v.isdigit() and int(v) >= 1)]
if bad_n:
    raise ValueError(f"sample sheet: min_ind_per_pop must be a positive integer (species {bad_n})")
bad = [s for s in SAMPLES.index if not re.fullmatch(r"[A-Za-z0-9_.-]+", s)]
if bad:
    raise ValueError(f"species names may only contain letters, digits, '_', '.', '-': {bad}")
if SAMPLES.index.duplicated().any():
    raise ValueError("duplicate species names in sample sheet")
# each species has its own popmap and latlongs files, so adding a species never
# means editing a file another species reads
for col in ("popmap", "latlongs"):
    paths = SAMPLES[col].map(os.path.abspath)
    shared = paths[paths.duplicated(keep=False)]
    if len(shared):
        raise ValueError(
            f"sample sheet: each species needs its own {col} file; shared by "
            + "; ".join(f"{p}: {', '.join(shared.index[shared == p])}" for p in shared.unique())
        )

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


_REFIT = config["gf"].get("refit_significant") or {}
REFIT_ON = bool(_REFIT.get("enabled", False))
REFIT_DIR = str(_REFIT.get("subdir", "refit_sigvars"))
if not re.fullmatch(r"[A-Za-z0-9_.-]+", REFIT_DIR):
    raise ValueError(f"gf.refit_significant.subdir must be a plain directory name, got {REFIT_DIR!r}")

_SP_RE = "|".join(re.escape(s) for s in SPECIES)

wildcard_constraints:
    species=_SP_RE,
    # a fitted model: "<species>" (all predictors) or "<species>/<REFIT_DIR>" (significant predictors)
    model=f"({_SP_RE})(/{re.escape(REFIT_DIR)})?",
    batch=r"\d+",


def min_ind_per_pop(wc):
    """per-species value from the sample sheet, else genotypes.min_ind_per_pop"""
    v = SAMPLES.loc[wc.species, "min_ind_per_pop"]
    return int(v) if v else int(config["genotypes"]["min_ind_per_pop"])


def sample_col(col):
    return lambda wc: SAMPLES.loc[wc.species, col]


def sp_dir(sub=""):
    return f"{OUT}/{{species}}/{sub}".rstrip("/")


def m_dir(sub=""):
    """path inside a model directory (wildcard {model})"""
    return f"{OUT}/{{model}}/{sub}".rstrip("/")


def model_species(wc):
    return wc.model.split("/")[0]


def is_refit(wc):
    return "/" in wc.model


def sp_file(sub):
    """input function: species-level file for a model"""
    return lambda wc: f"{OUT}/{model_species(wc)}/{sub}"


def model_predictors(wc):
    sp = model_species(wc)
    return f"{OUT}/{sp}/{REFIT_DIR}/predictors.txt" if is_refit(wc) else f"{OUT}/{sp}/predictors/selected.txt"


def model_label(wc):
    sp = model_species(wc)
    return f"{sp} (significant predictors only)" if is_refit(wc) else sp


def models():
    return SPECIES + ([f"{sp}/{REFIT_DIR}" for sp in SPECIES] if REFIT_ON else [])


def qc_html(model):
    """QC report name: qc.html for a full model, qc_refit.html for a refit"""
    return f"{OUT}/{model}/{'qc_refit' if '/' in model else 'qc'}.html"


def all_outputs():
    plots = [
        "map_turnover_rgb",
        "biplot_turnover_pca",
        "importance_r2weighted",
        "importance_accuracy",
        "random_vs_real",
    ]
    out = [f"{OUT}/{sp}/predictors/correlation.png" for sp in SPECIES]
    for m in models():
        for p in plots:
            out += [f"{OUT}/{m}/plots/{p}.pdf", f"{OUT}/{m}/plots/{p}.png"]
        out += [
            f"{OUT}/{m}/plots/gf_diagnostics.pdf",
            qc_html(m),
            f"{OUT}/{m}/gf/random_pvalues.tsv",
            f"{OUT}/{m}/gf/snp_pvalues.tsv",
        ]
    return out
