# Gradient forest fitting, randomization tests and spatial prediction.

GF = config["gf"]


localrules:
    install_gf,


rule install_gf:
    """One-time install of extendedForest + gradientForest (R-Forge) into the R conda env.
    Runs on the submit host because compute nodes may lack internet access."""
    output:
        f"{OUT}/.gf_installed",
    log:
        f"{OUT}/logs/install_gf.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/install_gf.R > {log} 2>&1 && touch {output}"


GF_ARGS = (
    f"--maxLevel={GF['maxLevel']} "
    f"--corr_threshold={GF['corr_threshold']} --nbin={GF['nbin']} --seed={GF['seed']}"
)
_TUNE = GF.get("ntree_tuning", {})
NTREE_ARGS = f"--ntree={GF.get('ntree', 'auto')}" + (
    (
        f" --tune_start={_TUNE.get('start', 250)} --tune_factor={_TUNE.get('factor', 2)}"
        f" --tune_max={_TUNE.get('max', 4000)} --min_snp_rsq_cor={_TUNE.get('min_snp_rsq_cor', 0.98)}"
        f" --max_importance_change={_TUNE.get('max_importance_change', 0.05)}"
        f" --min_turnover_cor={_TUNE.get('min_turnover_cor', 0.99)} --n_cells={_TUNE.get('n_cells', 1500)}"
    )
    if str(GF.get("ntree", "auto")) == "auto"
    else ""
)


rule gf_fit:
    input:
        gf_ok=rules.install_gf.output,
        sites=sp_file("env_sites.tsv"),
        freq=sp_file("genotypes/popfreq.tsv"),
        preds=model_predictors,
        tif=rules.prep_env_stack.output.tif,
    output:
        model=m_dir("gf/gf_model.rds"),
        importance=m_dir("gf/importance.tsv"),
        rsq=m_dir("gf/snp_rsq.tsv"),
        summary=m_dir("gf/summary.tsv"),
        tuning=m_dir("gf/ntree_tuning.tsv"),
        ntree=m_dir("gf/ntree_selected.tsv"),
    params:
        args=GF_ARGS,
        ntree=NTREE_ARGS,
    resources:
        mem_mb=24000,
        runtime=2880,  # adaptive ntree may fit up to ~2x tune_max trees in total
    log:
        f"{OUT}/logs/{{model}}/gf_fit.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_fit.R --env_sites={input.sites} --popfreq={input.freq} "
        "--predictors={input.preds} --env_tif={input.tif} {params.args} {params.ntree} "
        "--out_model={output.model} --out_importance={output.importance} "
        "--out_rsq={output.rsq} --out_summary={output.summary} "
        "--out_tuning={output.tuning} --out_ntree={output.ntree} > {log} 2>&1"


rule gf_randomize:
    input:
        gf_ok=rules.install_gf.output,
        sites=sp_file("env_sites.tsv"),
        freq=sp_file("genotypes/popfreq.tsv"),
        preds=model_predictors,
        ntree=m_dir("gf/ntree_selected.tsv"),
    output:
        stats=m_dir("gf/random/batch{batch}.tsv"),
        snp_rsq=m_dir("gf/random/batch{batch}_snp_rsq.tsv.gz"),
    params:
        args=GF_ARGS,
        first=lambda wc: batch_reps(wc)[0],
        last=lambda wc: batch_reps(wc)[1],
    resources:
        mem_mb=16000,
        runtime=random_runtime,
    log:
        f"{OUT}/logs/{{model}}/gf_random/batch{{batch}}.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_randomize.R --env_sites={input.sites} --popfreq={input.freq} "
        "--predictors={input.preds} --ntree_file={input.ntree} {params.args} "
        "--rep_start={params.first} --rep_end={params.last} "
        "--out={output.stats} --out_snp_rsq={output.snp_rsq} > {log} 2>&1"


rule gf_random_summary:
    input:
        reps=expand(m_dir("gf/random/batch{batch}.tsv"), batch=range(1, N_BATCHES + 1), allow_missing=True),
        snp_rsq=expand(
            m_dir("gf/random/batch{batch}_snp_rsq.tsv.gz"), batch=range(1, N_BATCHES + 1), allow_missing=True
        ),
        summary=m_dir("gf/summary.tsv"),
        importance=m_dir("gf/importance.tsv"),
    output:
        null=m_dir("gf/random_null.tsv"),
        pvals=m_dir("gf/random_pvalues.tsv"),
    log:
        f"{OUT}/logs/{{model}}/gf_random_summary.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_random_summary.R --summary={input.summary} "
        "--importance={input.importance} --out_null={output.null} --out_pvals={output.pvals} "
        "{input.reps} > {log} 2>&1"


rule gf_snp_pvalues:
    """per-SNP empirical p-values from the randomized fits, plus Benjamini-Hochberg FDR"""
    input:
        rsq=m_dir("gf/snp_rsq.tsv"),
        freq=sp_file("genotypes/popfreq.tsv"),
        null=m_dir("gf/random_null.tsv"),
        snp_rsq=expand(
            m_dir("gf/random/batch{batch}_snp_rsq.tsv.gz"), batch=range(1, N_BATCHES + 1), allow_missing=True
        ),
    output:
        m_dir("gf/snp_pvalues.tsv"),
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        f"{OUT}/logs/{{model}}/gf_snp_pvalues.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_snp_pvalues.R --rsq={input.rsq} --popfreq={input.freq} "
        "--null={input.null} --out={output} {input.snp_rsq} > {log} 2>&1"


rule gf_predict:
    input:
        gf_ok=rules.install_gf.output,
        model=m_dir("gf/gf_model.rds"),
        tif=rules.prep_env_stack.output.tif,
        sites=sp_file("env_sites.tsv"),
    output:
        turnover=m_dir("gf/turnover_stack.tif"),
        rgb=m_dir("gf/turnover_rgb.tif"),
        pca=m_dir("gf/pca.rds"),
        grid=m_dir("gf/grid_pcs.rds"),
        site_pcs=m_dir("gf/site_pcs.tsv"),
        loadings=m_dir("gf/pca_loadings.tsv"),
    params:
        pca_vars=GF.get("pca_vars", "all"),
    resources:
        mem_mb=48000,
        runtime=120,
    log:
        f"{OUT}/logs/{{model}}/gf_predict.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_predict.R --model={input.model} --env_tif={input.tif} "
        "--env_sites={input.sites} --pca_vars={params.pca_vars} "
        "--out_turnover={output.turnover} --out_rgb={output.rgb} --out_pca={output.pca} "
        "--out_grid={output.grid} --out_site_pcs={output.site_pcs} "
        "--out_loadings={output.loadings} > {log} 2>&1"


rule select_sig_predictors:
    """Predictors whose observed importance beats the randomized fits (full model)."""
    input:
        pvals=sp_dir("gf/random_pvalues.tsv"),
    output:
        selected=f"{OUT}/{{species}}/{REFIT_DIR}/predictors.txt",
        report=f"{OUT}/{{species}}/{REFIT_DIR}/predictor_selection.tsv",
    params:
        alpha=_REFIT.get("alpha", 0.05),
        metrics=",".join(_REFIT.get("metrics", ["importance", "accuracy"])),
        min_n=_REFIT.get("min_predictors", 2),
    log:
        f"{OUT}/logs/{{species}}/select_sig_predictors.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/select_sig_predictors.py --pvals {input.pvals} --alpha {params.alpha} "
        "--metrics {params.metrics} --min-predictors {params.min_n} "
        "--out-list {output.selected} --out-report {output.report} > {log} 2>&1"
