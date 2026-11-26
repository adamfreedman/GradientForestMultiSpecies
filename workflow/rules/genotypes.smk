# Input validation, VCF -> population allele frequencies, and a genotype PCA.


rule validate_inputs:
    input:
        vcf=sample_col("vcf"),
        popmap=sample_col("popmap"),
        latlongs=sample_col("latlongs"),
        tif=rules.prep_env_stack.output.tif,
    output:
        report=sp_dir("qc/validation.txt"),
        sample_pop=sp_dir("qc/sample_pop.tsv"),
        unmapped=sp_dir("qc/unmapped_samples.tsv"),
    params:
        regex=lambda wc: SAMPLES.loc[wc.species, "id_regex"],
        snap=config["env"].get("max_snap_km", 0),
    log:
        f"{OUT}/logs/{{species}}/validate_inputs.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/validate_inputs.py --vcf {input.vcf} --popmap {input.popmap} "
        "--latlongs {input.latlongs} --env-tif {input.tif} --id-regex '{params.regex}' "
        "--max-snap-km {params.snap} --out-report {output.report} "
        "--out-sample-pop {output.sample_pop} --out-unmapped {output.unmapped} > {log} 2>&1"


rule vcf_to_popfreq:
    input:
        vcf=sample_col("vcf"),
        sample_pop=sp_dir("qc/sample_pop.tsv"),
        sites=sp_dir("env_sites.tsv"),
    output:
        freq=sp_dir("genotypes/popfreq.tsv"),
        summary=sp_dir("genotypes/snp_filter_summary.tsv"),
        pops=sp_dir("genotypes/pop_counts.tsv"),
        snp_info=sp_dir("genotypes/snp_info.tsv"),
    params:
        min_ind=min_ind_per_pop,
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        f"{OUT}/logs/{{species}}/vcf_to_popfreq.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/vcf_to_popfreq.py --vcf {input.vcf} --sample-pop {input.sample_pop} "
        "--env-sites {input.sites} --min-ind-per-pop {params.min_ind} "
        "--out {output.freq} --out-summary {output.summary} --out-pops {output.pops} "
        "--out-snp-info {output.snp_info} > {log} 2>&1"


rule genotype_pca:
    """PCA of individual genotypes: population structure and sample QC (not used by GF)."""
    input:
        vcf=sample_col("vcf"),
        sample_pop=sp_dir("qc/sample_pop.tsv"),
        sites=sp_dir("env_sites.tsv"),
        pops=sp_dir("genotypes/pop_counts.tsv"),
    output:
        variance=sp_dir("pca/variance_explained.tsv"),
        scores=sp_dir("pca/sample_scores.tsv"),
        info=sp_dir("pca/pca_info.tsv"),
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        f"{OUT}/logs/{{species}}/genotype_pca.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/genotype_pca.py --vcf {input.vcf} --sample-pop {input.sample_pop} "
        "--env-sites {input.sites} --pop-counts {input.pops} --out-variance {output.variance} "
        "--out-scores {output.scores} --out-info {output.info} > {log} 2>&1"


rule plot_genotype_pca:
    input:
        scores=sp_dir("pca/sample_scores.tsv"),
        variance=sp_dir("pca/variance_explained.tsv"),
        info=sp_dir("pca/pca_info.tsv"),
    output:
        pdf=sp_dir("pca/pca_pc1_pc2.pdf"),
        png=sp_dir("pca/pca_pc1_pc2.png"),
    params:
        prefix=lambda wc, output: output.pdf[: -len(".pdf")],
    log:
        f"{OUT}/logs/{{species}}/plot_genotype_pca.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/plot_genotype_pca.R --species='{wildcards.species}' --scores={input.scores} "
        "--variance={input.variance} --info={input.info} --out_prefix={params.prefix} > {log} 2>&1"
