# ADMIXTURE ancestry proportions over a range of K, cross-validation to choose K,
# and pong to check that replicate runs at each K agree (population structure; not
# used by gradient forest). Settings: admixture in config.yaml.


rule admixture_input:
    """VCF -> PLINK bed for ADMIXTURE: the genotype PCA's individuals (every sample mapped
    to a population with coordinates) and biallelic SNPs variable among them. ADMIXTURE
    only accepts integer chromosome codes; it does not use them, so all are set to 1."""
    input:
        vcf=sample_col("vcf"),
        sample_pop=sp_dir("qc/sample_pop.tsv"),
    output:
        bed=sp_dir("admixture/input/geno.bed"),
        bim=sp_dir("admixture/input/geno.bim"),
        fam=sp_dir("admixture/input/geno.fam"),
        keep=sp_dir("admixture/input/keep.txt"),
    params:
        prefix=sp_dir("admixture/input/geno"),
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        f"{OUT}/logs/{{species}}/admixture_input.log",
    conda:
        "../envs/admixture.yaml"
    shell:
        "awk 'NR > 1 {{print $1}}' {input.sample_pop} > {output.keep} && "
        "plink2 --vcf {input.vcf} --keep {output.keep} --snps-only just-acgt --max-alleles 2 --mac 1 "
        "--vcf-half-call missing --allow-extra-chr --set-all-var-ids '@:#' --make-bed "
        "--threads {threads} --memory 6000 --out {params.prefix} > {log} 2>&1 && "
        "awk 'BEGIN {{OFS = \"\\t\"}} {{$1 = 1; print}}' {output.bim} > {output.bim}.tmp && "
        "mv {output.bim}.tmp {output.bim}"


rule admixture_run:
    """one ADMIXTURE run: K clusters, replicate rep (seed = admixture.seed + rep - 1),
    with cross-validation"""
    input:
        bed=rules.admixture_input.output.bed,
        bim=rules.admixture_input.output.bim,
        fam=rules.admixture_input.output.fam,
    output:
        q=sp_dir("admixture/runs/K{k}/rep{rep}/geno.{k}.Q"),
        p=sp_dir("admixture/runs/K{k}/rep{rep}/geno.{k}.P"),
        log=sp_dir("admixture/runs/K{k}/rep{rep}/admixture.log"),
    params:
        dir=sp_dir("admixture/runs/K{k}/rep{rep}"),
        abs_in=lambda wc, input: " ".join(os.path.abspath(f) for f in (input.bed, input.bim, input.fam)),
        seed=lambda wc: ADMIX_SEED + int(wc.rep) - 1,
        cv=ADMIX_CV,
    threads: ADMIX_THREADS
    resources:
        mem_mb=4000,
        # ~3.5 min for 208 individuals, 56k SNPs, K = 10, 10-fold CV, 4 threads
        runtime=120,
    conda:
        "../envs/admixture.yaml"
    shell:
        "cd {params.dir} && ln -sf {params.abs_in} . && "
        "admixture --cv={params.cv} -s {params.seed} -j{threads} geno.bed {wildcards.k} > admixture.log 2>&1 && "
        "rm -f geno.bed geno.bim geno.fam"


def admixture_runs(wc, kind="q"):
    name = {"q": "geno.{k}.Q", "log": "admixture.log"}[kind]
    return [
        f"{OUT}/{wc.species}/admixture/runs/K{k}/rep{r}/" + name.format(k=k)
        for k in ADMIX_K
        for r in range(1, ADMIX_REPS + 1)
    ]


rule admixture_pong:
    """pong on every run: groups the replicate runs at each K into modes"""
    input:
        q=admixture_runs,
        fam=rules.admixture_input.output.fam,
        sample_pop=sp_dir("qc/sample_pop.tsv"),
    output:
        summary=sp_dir("admixture/pong/output/result_summary.txt"),
        filemap=sp_dir("admixture/pong/filemap.txt"),
    params:
        k=",".join(map(str, ADMIX_K)),
        reps=ADMIX_REPS,
        sim=ADMIX_SIM,
    log:
        f"{OUT}/logs/{{species}}/admixture_pong.log",
    conda:
        "../envs/admixture.yaml"
    shell:
        "python {SCRIPTS}/admixture_pong.py --runs-dir {OUT}/{wildcards.species}/admixture/runs "
        "--k {params.k} --reps {params.reps} --fam {input.fam} --sample-pop {input.sample_pop} "
        "--sim-threshold {params.sim} --filemap {output.filemap} "
        "--outdir $(dirname {output.summary}) > {log} 2>&1"


rule admixture_summary:
    """CV error per K, pong modes, choice of K, and ancestry tables with clusters
    aligned across K"""
    input:
        q=admixture_runs,
        logs=lambda wc: admixture_runs(wc, "log"),
        fam=rules.admixture_input.output.fam,
        sample_pop=sp_dir("qc/sample_pop.tsv"),
        pong=rules.admixture_pong.output.summary,
    output:
        cv=sp_dir("admixture/cv_error.tsv"),
        k_summary=sp_dir("admixture/k_summary.tsv"),
        selected=sp_dir("admixture/selected_k.tsv"),
        modes=sp_dir("admixture/selected_k_modes.tsv"),
        ancestry=expand(sp_dir("admixture/ancestry_K{k}.tsv"), k=ADMIX_K, allow_missing=True),
    params:
        k=",".join(map(str, ADMIX_K)),
        reps=ADMIX_REPS,
        seed=ADMIX_SEED,
        frac=ADMIX_STABLE_FRAC,
        outdir=sp_dir("admixture"),
    log:
        f"{OUT}/logs/{{species}}/admixture_summary.log",
    conda:
        "../envs/admixture.yaml"
    shell:
        "python {SCRIPTS}/admixture_summary.py --runs-dir {OUT}/{wildcards.species}/admixture/runs "
        "--k {params.k} --reps {params.reps} --seed {params.seed} --fam {input.fam} "
        "--sample-pop {input.sample_pop} --pong-summary {input.pong} --stable-min-frac {params.frac} "
        "--outdir {params.outdir} > {log} 2>&1"


rule plot_admixture:
    input:
        cv=rules.admixture_summary.output.cv,
        k_summary=rules.admixture_summary.output.k_summary,
        selected=rules.admixture_summary.output.selected,
        modes=rules.admixture_summary.output.modes,
        ancestry=rules.admixture_summary.output.ancestry,
    output:
        expand(sp_dir("admixture/{plot}.{ext}"), plot=ADMIX_PLOTS, ext=["pdf", "png"], allow_missing=True),
    params:
        outdir=sp_dir("admixture"),
    log:
        f"{OUT}/logs/{{species}}/plot_admixture.log",
    conda:
        "../envs/r_gf.yaml"
    shell:
        "Rscript {SCRIPTS}/plot_admixture.R --species='{wildcards.species}' --admixdir={params.outdir} "
        "> {log} 2>&1"
