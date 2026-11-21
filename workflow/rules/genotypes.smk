# Input validation and VCF -> population allele frequencies.


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
    params:
        min_ind=config["genotypes"]["min_ind_per_pop"],
        max_miss=config["genotypes"]["max_missing_per_pop"],
        maf=config["genotypes"]["min_maf"],
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
        "--max-missing-per-pop {params.max_miss} --min-maf {params.maf} "
        "--out {output.freq} --out-summary {output.summary} --out-pops {output.pops} > {log} 2>&1"
