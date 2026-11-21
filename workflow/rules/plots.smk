# Single-species figures.


rule plot_gf:
    input:
        gf_ok=rules.install_gf.output,
        model=m_dir("gf/gf_model.rds"),
        importance=m_dir("gf/importance.tsv"),
        summary=m_dir("gf/summary.tsv"),
        null=m_dir("gf/random_null.tsv"),
        pvals=m_dir("gf/random_pvalues.tsv"),
        grid=m_dir("gf/grid_pcs.rds"),
        pca=m_dir("gf/pca.rds"),
        site_pcs=m_dir("gf/site_pcs.tsv"),
        area=rules.study_area.output.geojson,
        boundaries=rules.get_boundaries.output,
    output:
        expand(
            m_dir("plots/{plot}.{ext}"),
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
        diag=m_dir("plots/gf_diagnostics.pdf"),
    params:
        outdir=m_dir("plots"),
        label=model_label,
        maxpts=config["gf"].get("biplot_max_points", 50000),
    resources:
        mem_mb=32000,
        runtime=60,
    log:
        f"{OUT}/logs/{{model}}/plot_gf.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/plot_gf.R --species='{params.label}' --model={input.model} "
        "--importance={input.importance} --summary={input.summary} --null={input.null} "
        "--pvals={input.pvals} --grid={input.grid} --pca={input.pca} "
        "--site_pcs={input.site_pcs} --study_area={input.area} --boundaries={input.boundaries} "
        "--biplot_max_points={params.maxpts} --outdir={params.outdir} > {log} 2>&1"


rule qc_report:
    input:
        plots=rules.plot_gf.output,
        validation=sp_file("qc/validation.txt"),
        snp_filter=sp_file("genotypes/snp_filter_summary.tsv"),
        pop_counts=sp_file("genotypes/pop_counts.tsv"),
        sites=sp_file("env_sites.tsv"),
        predictors=lambda wc: (
            f"{OUT}/{wc.model}/predictor_selection.tsv" if is_refit(wc) else f"{OUT}/{wc.model}/predictors/report.tsv"
        ),
        heatmap=sp_file("predictors/correlation.png"),
        tuning=m_dir("gf/ntree_tuning.tsv"),
        ntree=m_dir("gf/ntree_selected.tsv"),
        summary=m_dir("gf/summary.tsv"),
        pvals=m_dir("gf/random_pvalues.tsv"),
    output:
        html=m_dir("qc.html"),
        tuning_png=m_dir("plots/ntree_tuning.png"),
    params:
        spdir=m_dir(),
        label=model_label,
        refit=lambda wc: "--refit=1" if is_refit(wc) else "",
        thr=lambda wc: (
            f"--min_snp_rsq_cor={_TUNE.get('min_snp_rsq_cor', 0.98)} "
            f"--max_importance_change={_TUNE.get('max_importance_change', 0.05)} "
            f"--min_turnover_cor={_TUNE.get('min_turnover_cor', 0.99)}"
        ),
    log:
        f"{OUT}/logs/{{model}}/qc_report.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/qc_report.R --species='{params.label}' --spdir={params.spdir} "
        "--validation={input.validation} --snp_filter={input.snp_filter} --pop_counts={input.pop_counts} "
        "--env_sites={input.sites} --predictors={input.predictors} --tuning={input.tuning} "
        "--ntree={input.ntree} --summary={input.summary} --pvals={input.pvals} {params.thr} "
        "--heatmap={input.heatmap} {params.refit} "
        "--out={output.html} --out_tuning_png={output.tuning_png} > {log} 2>&1"
