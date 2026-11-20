# Study area, environmental raster stack and per-site environment.


def study_area_args():
    clip = config["env"].get("clip") or ""
    bbox = config["env"].get("bbox") or []
    if clip:
        return f"--vector '{clip}'"
    if len(bbox) == 4:
        return "--extent " + " ".join(str(b) for b in bbox)
    raise ValueError("config env: set either 'clip' (vector file) or 'bbox' [xmin, ymin, xmax, ymax]")


rule study_area:
    output:
        geojson=f"{ENV_DIR}/study_area.geojson",
        bbox=f"{ENV_DIR}/study_area_bbox.tsv",
    params:
        src=study_area_args(),
    log:
        f"{OUT}/logs/study_area.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/bbox_to_vector.py {params.src} "
        "--geojson {output.geojson} --bbox {output.bbox} > {log} 2>&1"


localrules:
    get_boundaries,


rule get_boundaries:
    """Country boundaries clipped to the study area (runs on the submit host: needs internet)."""
    input:
        geojson=rules.study_area.output.geojson,
    output:
        f"{ENV_DIR}/boundaries.geojson",
    params:
        source=lambda wc: config["env"].get("boundaries") or "none",
        cache=f"{ENV_DIR}/downloads",
    log:
        f"{OUT}/logs/get_boundaries.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/get_boundaries.py --source '{params.source}' --clip {input.geojson} "
        "--cache-dir {params.cache} --out {output} > {log} 2>&1"


rule prep_env_stack:
    input:
        geojson=rules.study_area.output.geojson,
    output:
        tif=f"{ENV_DIR}/env_stack.tif",
        layers=f"{ENV_DIR}/layers.txt",
    params:
        dir=config["env"]["dir"],
        pattern=config["env"]["pattern"],
        # lambdas stop Snakemake from treating "{0}" etc. as wildcards
        regex=lambda wc: config["env"]["name_regex"],
        template=lambda wc: config["env"].get("name_template", "bio{0}"),
        agg=config["env"].get("grid_aggregate", 1),
    resources:
        mem_mb=24000,
    log:
        f"{OUT}/logs/prep_env_stack.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/prep_env_stack.py --env-dir '{params.dir}' --pattern '{params.pattern}' "
        "--name-regex '{params.regex}' --name-template '{params.template}' "
        "--clip {input.geojson} --aggregate {params.agg} "
        "--out-tif {output.tif} --out-layers {output.layers} > {log} 2>&1"


rule extract_env:
    input:
        latlongs=sample_col("latlongs"),
        sample_pop=sp_dir("qc/sample_pop.tsv"),
        tif=rules.prep_env_stack.output.tif,
    output:
        sp_dir("env_sites.tsv"),
    params:
        buffer=config["env"].get("site_buffer_km", 0),
        snap=config["env"].get("max_snap_km", 0),
    log:
        f"{OUT}/logs/{{species}}/extract_env.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/extract_env.py --latlongs {input.latlongs} --sample-pop {input.sample_pop} "
        "--env-tif {input.tif} --buffer-km {params.buffer} --max-snap-km {params.snap} "
        "--out {output} > {log} 2>&1"


def predictor_include():
    inc = config["predictors"].get("include", "all")
    return "all" if inc in (None, "all") else ",".join(inc)


rule select_predictors:
    input:
        sites=sp_dir("env_sites.tsv"),
        layers=rules.prep_env_stack.output.layers,
    output:
        selected=sp_dir("predictors/selected.txt"),
        report=sp_dir("predictors/report.tsv"),
        heatmap=sp_dir("predictors/correlation.png"),
    params:
        include=predictor_include(),
        thr=config["predictors"].get("corr_threshold") or "none",
        keep=",".join(config["predictors"].get("keep_priority") or []) or "none",
    log:
        f"{OUT}/logs/{{species}}/select_predictors.log",
    conda:
        "../envs/py_geo.yaml"
    shell:
        "python {SCRIPTS}/select_predictors.py --env-sites {input.sites} --layers {input.layers} "
        "--include {params.include} --threshold {params.thr} --keep-priority {params.keep} "
        "--out-list {output.selected} --out-report {output.report} --out-heatmap {output.heatmap} "
        "> {log} 2>&1"
