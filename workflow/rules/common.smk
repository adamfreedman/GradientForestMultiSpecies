import math
import os
import re

import pandas as pd

SCRIPTS = os.path.join(workflow.basedir, "scripts")
def _flag(x):
    return str(x).strip().lower() in ("true", "1", "yes")


# use_all_predictors (run_slurm.sh --use-all-predictors, or --config use_all_predictors=true):
# no predictor selection; every candidate variable goes into the model. Results then go to
# <outdir>_all_predictors so they never overwrite a run with selection. The environmental
# layers and the gradientForest install are shared (always under <outdir>).
USE_ALL_PREDICTORS = _flag(config.get("use_all_predictors", False))
# importance_floor (run_slurm.sh --importance-floor 0.01, or --config importance_floor=0.01):
# drop selected variables whose share of the summed screening importance is below this
# fraction. Default: no floor. Results then go to <outdir>_importance_floor<value>.
_floor = config.get("importance_floor")
IMPORTANCE_FLOOR = None if _floor in (None, "", "none", "null") else float(_floor)
if IMPORTANCE_FLOOR is not None:
    if not 0 <= IMPORTANCE_FLOOR < 1:
        raise ValueError(f"importance_floor must be a fraction in [0, 1), got {_floor!r}")
    if USE_ALL_PREDICTORS:
        raise ValueError("importance_floor applies to predictor selection; it cannot be combined with use_all_predictors")
BASE_OUT = config.get("outdir", "results")
if USE_ALL_PREDICTORS:
    OUT = f"{BASE_OUT}_all_predictors"
elif IMPORTANCE_FLOOR is not None:
    OUT = f"{BASE_OUT}_importance_floor{IMPORTANCE_FLOOR:g}"
else:
    OUT = BASE_OUT
ENV_DIR = f"{BASE_OUT}/env"

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


def snp_scale(freq):
    """number of SNP columns in popfreq.tsv / 14,000 (the size the base
    resources were set for), at least 1"""
    with open(freq) as fh:
        n = fh.readline().count("\t")
    return max(1.0, n / 14000)


_PRED = config["predictors"]
SELECT_THRESHOLD = "none" if USE_ALL_PREDICTORS else (_PRED.get("corr_threshold") or "none")
# leave-one-population-out predictor selection (diagnostic table in qc.html); only
# meaningful when variables are selected
STABILITY_ON = _flag(_PRED.get("stability", True)) and SELECT_THRESHOLD != "none"
SELECT_RANK_BY = _PRED.get("rank_by", "importance")
if SELECT_RANK_BY not in ("importance", "accuracy", "correlation"):
    raise ValueError(f"predictors.rank_by must be importance, accuracy or correlation, got {SELECT_RANK_BY!r}")
SCREEN_NTREE = int(_PRED.get("screen_ntree", 500))
SELECT_ARGS = (
    f"--rank_by={SELECT_RANK_BY} --screen_ntree={SCREEN_NTREE} "
    f"--keep_priority={','.join(_PRED.get('keep_priority') or []) or 'none'}"
)
if IMPORTANCE_FLOOR is not None:
    if SELECT_THRESHOLD == "none" or SELECT_RANK_BY == "correlation":
        raise ValueError("importance_floor needs predictors.corr_threshold set and rank_by importance or accuracy")
    # only added when set, so runs without a floor keep their params (no reruns)
    SELECT_ARGS += f" --importance_floor={IMPORTANCE_FLOOR:g}"


def random_runtime(wc, input):
    """SLURM minutes for one randomization batch. Measured ~1.1 min per 100 trees
    per fit for 56k SNPs, 9 predictors and 14 sites; allow 1 min per 100 trees
    per 10 predictors per 14k SNPs (~3x margin)."""
    ntree = int(pd.read_csv(input.ntree, sep="\t")["ntree"].iloc[0])
    n_pred = sum(1 for l in open(input.preds) if l.strip())
    per_rep = ntree / 100 * max(1.0, n_pred / 10) * snp_scale(input.freq)
    return 30 + RANDOM_BATCH * math.ceil(per_rep)


def n_populations(sp):
    """populations in the species' latlongs file: an upper bound on the populations
    with genotype data, so the number of leave-one-out jobs is known before the
    genotypes are processed (jobs past the real number write empty tables)"""
    return len(pd.read_csv(SAMPLES.loc[sp, "latlongs"], dtype=str))


def screen_runtime(wc, input):
    """SLURM minutes for one screening fit with every candidate: measured ~18 min at
    500 trees, 19 variables, 56k SNPs, 13 sites; ~4x margin"""
    n_pred = sum(1 for l in open(input.candidates) if l.strip())
    return 30 + math.ceil(SCREEN_NTREE / 100 * max(1.0, n_pred / 10) * snp_scale(input.freq))


def random_mem(wc, input):
    return math.ceil(16000 * snp_scale(input.freq))


def fit_mem(wc, input):
    return math.ceil(24000 * snp_scale(input.freq))


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
    out += [f"{OUT}/{sp}/pca/pca_pc1_pc2.{ext}" for sp in SPECIES for ext in ("pdf", "png")]
    if STABILITY_ON:
        out += [f"{OUT}/{sp}/predictors/stability_variables.tsv" for sp in SPECIES]
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
