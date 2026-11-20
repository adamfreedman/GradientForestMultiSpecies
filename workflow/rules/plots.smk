# Single-species figures.


rule plot_gf:
    input:
        gf_ok=rules.install_gf.output,
        model=sp_dir("gf/gf_model.rds"),
        importance=sp_dir("gf/importance.tsv"),
        summary=sp_dir("gf/summary.tsv"),
        null=sp_dir("gf/random_null.tsv"),
        pvals=sp_dir("gf/random_pvalues.tsv"),
        grid=sp_dir("gf/grid_pcs.rds"),
        pca=sp_dir("gf/pca.rds"),
        site_pcs=sp_dir("gf/site_pcs.tsv"),
        area=rules.study_area.output.geojson,
        boundaries=rules.get_boundaries.output,
    output:
        expand(
            sp_dir("plots/{plot}.{ext}"),
            plot=[
                "map_turnover_rgb",
                "biplot_turnover_pca",
                "importance_r2weighted",
                "importance_accuracy",
                "random_vs_real",
            ],
            ext=["pdf", "png"],
            allow_missing=True,
        ),
        diag=sp_dir("plots/gf_diagnostics.pdf"),
    params:
        outdir=sp_dir("plots"),
        maxpts=config["gf"].get("biplot_max_points", 50000),
    resources:
        mem_mb=32000,
        runtime=60,
    log:
        f"{OUT}/logs/{{species}}/plot_gf.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/plot_gf.R --species={wildcards.species} --model={input.model} "
        "--importance={input.importance} --summary={input.summary} --null={input.null} "
        "--pvals={input.pvals} --grid={input.grid} --pca={input.pca} "
        "--site_pcs={input.site_pcs} --study_area={input.area} --boundaries={input.boundaries} "
        "--biplot_max_points={params.maxpts} --outdir={params.outdir} > {log} 2>&1"


rule qc_report:
    input:
        plots=rules.plot_gf.output,
        validation=sp_dir("qc/validation.txt"),
        snp_filter=sp_dir("genotypes/snp_filter_summary.tsv"),
        pop_counts=sp_dir("genotypes/pop_counts.tsv"),
        sites=sp_dir("env_sites.tsv"),
        predictors=sp_dir("predictors/report.tsv"),
        heatmap=sp_dir("predictors/correlation.png"),
        tuning=sp_dir("gf/ntree_tuning.tsv"),
        ntree=sp_dir("gf/ntree_selected.tsv"),
        summary=sp_dir("gf/summary.tsv"),
        pvals=sp_dir("gf/random_pvalues.tsv"),
    output:
        html=sp_dir("qc.html"),
        tuning_png=sp_dir("plots/ntree_tuning.png"),
    params:
        spdir=sp_dir(),
        thr=lambda wc: (
            f"--min_snp_rsq_cor={_TUNE.get('min_snp_rsq_cor', 0.98)} "
            f"--max_importance_change={_TUNE.get('max_importance_change', 0.05)} "
            f"--min_turnover_cor={_TUNE.get('min_turnover_cor', 0.99)}"
        ),
    log:
        f"{OUT}/logs/{{species}}/qc_report.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/qc_report.R --species={wildcards.species} --spdir={params.spdir} "
        "--validation={input.validation} --snp_filter={input.snp_filter} --pop_counts={input.pop_counts} "
        "--env_sites={input.sites} --predictors={input.predictors} --tuning={input.tuning} "
        "--ntree={input.ntree} --summary={input.summary} --pvals={input.pvals} {params.thr} "
        "--out={output.html} --out_tuning_png={output.tuning_png} > {log} 2>&1"
