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
        sites=sp_dir("env_sites.tsv"),
        freq=sp_dir("genotypes/popfreq.tsv"),
        preds=sp_dir("predictors/selected.txt"),
        tif=rules.prep_env_stack.output.tif,
    output:
        model=sp_dir("gf/gf_model.rds"),
        importance=sp_dir("gf/importance.tsv"),
        rsq=sp_dir("gf/snp_rsq.tsv"),
        summary=sp_dir("gf/summary.tsv"),
        tuning=sp_dir("gf/ntree_tuning.tsv"),
        ntree=sp_dir("gf/ntree_selected.tsv"),
    params:
        args=GF_ARGS,
        ntree=NTREE_ARGS,
    resources:
        mem_mb=24000,
        runtime=2880,  # adaptive ntree may fit up to ~2x tune_max trees in total
    log:
        f"{OUT}/logs/{{species}}/gf_fit.log",
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
        sites=sp_dir("env_sites.tsv"),
        freq=sp_dir("genotypes/popfreq.tsv"),
        preds=sp_dir("predictors/selected.txt"),
        ntree=sp_dir("gf/ntree_selected.tsv"),
    output:
        stats=sp_dir("gf/random/batch{batch}.tsv"),
        snp_rsq=sp_dir("gf/random/batch{batch}_snp_rsq.tsv.gz"),
    params:
        args=GF_ARGS,
        first=lambda wc: batch_reps(wc)[0],
        last=lambda wc: batch_reps(wc)[1],
    resources:
        mem_mb=16000,
        runtime=random_runtime,
    log:
        f"{OUT}/logs/{{species}}/gf_random/batch{{batch}}.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_randomize.R --env_sites={input.sites} --popfreq={input.freq} "
        "--predictors={input.preds} --ntree_file={input.ntree} {params.args} "
        "--rep_start={params.first} --rep_end={params.last} "
        "--out={output.stats} --out_snp_rsq={output.snp_rsq} > {log} 2>&1"


rule gf_random_summary:
    input:
        reps=expand(sp_dir("gf/random/batch{batch}.tsv"), batch=range(1, N_BATCHES + 1), allow_missing=True),
        snp_rsq=expand(
            sp_dir("gf/random/batch{batch}_snp_rsq.tsv.gz"), batch=range(1, N_BATCHES + 1), allow_missing=True
        ),
        summary=sp_dir("gf/summary.tsv"),
        importance=sp_dir("gf/importance.tsv"),
    output:
        null=sp_dir("gf/random_null.tsv"),
        pvals=sp_dir("gf/random_pvalues.tsv"),
    log:
        f"{OUT}/logs/{{species}}/gf_random_summary.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_random_summary.R --summary={input.summary} "
        "--importance={input.importance} --out_null={output.null} --out_pvals={output.pvals} "
        "{input.reps} > {log} 2>&1"


rule gf_predict:
    input:
        gf_ok=rules.install_gf.output,
        model=sp_dir("gf/gf_model.rds"),
        tif=rules.prep_env_stack.output.tif,
        sites=sp_dir("env_sites.tsv"),
    output:
        turnover=sp_dir("gf/turnover_stack.tif"),
        rgb=sp_dir("gf/turnover_rgb.tif"),
        pca=sp_dir("gf/pca.rds"),
        grid=sp_dir("gf/grid_pcs.rds"),
        site_pcs=sp_dir("gf/site_pcs.tsv"),
        loadings=sp_dir("gf/pca_loadings.tsv"),
    params:
        pca_vars=GF.get("pca_vars", "all"),
    resources:
        mem_mb=48000,
        runtime=120,
    log:
        f"{OUT}/logs/{{species}}/gf_predict.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/gf_predict.R --model={input.model} --env_tif={input.tif} "
        "--env_sites={input.sites} --pca_vars={params.pca_vars} "
        "--out_turnover={output.turnover} --out_rgb={output.rgb} --out_pca={output.pca} "
        "--out_grid={output.grid} --out_site_pcs={output.site_pcs} "
        "--out_loadings={output.loadings} > {log} 2>&1"
