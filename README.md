# GradientForestMultiSpecies

Snakemake workflow for gradient forest (GF) models of genomic turnover across
environmental gradients, for one or more species sampled in the same region.
It generalizes the single-species analysis in
[TrachylepisAffinisSpeciation/gradient_forest](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/tree/master/gradient_forest)
and writes per-species turnover rasters on a common grid, for later
cross-species congruence analyses.

## Inputs

`config/samples.tsv` (tab-separated, one row per species). Each species has
its own popmap and latlongs files, named with a species prefix, e.g.
`taffinis_id_to_pop.txt` and `taffinis_pop_latlongs.csv`. To add a species, add
its files and one row to the sample sheet; no shared file is edited. The
workflow stops if two species point at the same popmap or latlongs file (a shared multi-species VCF is fine).

| column | content |
|---|---|
| `species` | short name used in output paths (letters, digits, `_ . -`) |
| `vcf` | genotype VCF (plain or bgzipped) |
| `popmap` | whitespace-delimited, header `sampleid popid`: VCF sample → population |
| `latlongs` | CSV, header `popid,sitename,lat,long` |
| `id_regex` | optional column: regex deleted from VCF sample IDs before matching `popmap` (e.g. `^MCZA_`); omit it when the IDs already match |

`config/config.yaml` sets the environmental-layer directory, the study area
(any vector file, or a bounding box), predictor pruning, SNP filters and GF
parameters. Each option is commented in the file.

## Steps

| rule | language | output (under `results/`) |
|---|---|---|
| `study_area` | Python | `env/study_area.geojson`, `env/study_area_bbox.tsv` |
| `get_boundaries` | Python | `env/boundaries.geojson`: country borders for the map (Natural Earth 1:10m by default, downloaded once) |
| `prep_env_stack` | Python | `env/env_stack.tif` (clipped multiband stack, bands named `bio1..bio19`) |
| `validate_inputs` | Python | `{sp}/qc/validation.txt`, `{sp}/qc/sample_pop.tsv`, `{sp}/qc/unmapped_samples.tsv`. Individuals that do not map to a population with coordinates (in the VCF but not the popmap, in the popmap but not the VCF, or in a popid missing from latlongs) are excluded with a warning and listed in `unmapped_samples.tsv` and `qc.html`. Fails, naming the population(s) in the log, if a latlongs population has no individuals in the VCF; also fails on duplicate IDs or sites off the environmental grid |
| `extract_env` | Python | `{sp}/env_sites.tsv` (environment at each population) |
| `select_predictors` | Python | `{sp}/predictors/selected.txt`, `report.tsv`, `correlation.png` |
| `vcf_to_popfreq` | Python | `{sp}/genotypes/popfreq.tsv` (populations × SNPs, minor-allele frequency) |
| `install_gf` | R | one-time install of `extendedForest`/`gradientForest` from R-Forge |
| `gf_fit` | R | `{sp}/gf/gf_model.rds`, `importance.tsv`, `snp_rsq.tsv`, `summary.tsv` |
| `gf_randomize` | R | `{sp}/gf/random/batch{n}.tsv`: `n_random` fits (default 1000) with environment rows permuted and the same GF settings, run `random_batch_size` per SLURM job. `batch{n}_snp_rsq.tsv.gz` keeps each replicate's per-SNP R² for later per-SNP p-values |
| `gf_random_summary` | R | `{sp}/gf/random_null.tsv`, `random_pvalues.tsv` |
| `gf_snp_pvalues` | R | `{sp}/gf/snp_pvalues.tsv`: per-SNP observed R², empirical p = (1 + no. randomized fits with R² ≥ observed) / (n_random + 1) (SNPs with observed R² = 0 get p = 1), and `p_bh`, the Benjamini-Hochberg adjustment (`p.adjust(method = "BH")`) over all SNPs tested |
| `gf_predict` | R | `{sp}/gf/turnover_stack.tif` (transformed predictors), `turnover_rgb.tif`, PCA objects |
| `qc_report` | R | `{sp}/qc.html` (`{sp}/refit_sigvars/qc_refit.html` for the refit), one self-contained file with figures embedded: warnings (including ntree not converging), ntree selection, model and randomization summaries, figures, predictor and SNP filtering, populations |
| `plot_gf` | R | `{sp}/plots/`: turnover map, PCA biplot, R²-weighted and accuracy importance charts, observed-vs-random histograms, GF diagnostics |

### Converting the study-area shapefile without a GIS GUI

```bash
python workflow/scripts/bbox_to_vector.py --vector data/bounding_box/PIREbb.shp \
    --geojson study_area.geojson --bbox study_area_bbox.tsv
```

The GeoJSON is plain text and can be read by R (`terra::vect`), Python
(`geopandas`), or QGIS. The TSV has the `xmin ymin xmax ymax` extent.

## Running on Cannon

```bash
conda activate my_snakemake
snakemake -n                                  # dry run
snakemake --cores 1 --use-conda results/.gf_installed   # one-time R package install (login node)
snakemake --use-conda --cores 1 results/env/boundaries.geojson   # downloads country borders (login node)
sbatch run_slurm.sh                           # everything else: a controller job that submits each rule to SLURM
```

`run_slurm.sh` passes any extra arguments to Snakemake. The controller's log is `logs/snakemake_head_<jobid>.log`; per-rule logs are under `results/logs/`.

For a quick test, set `gf: n_random: 5` in the config.

## Notes on methods

* The responses are population allele frequencies (minor allele), with one row per population, as in the original analysis.
* `ntree: auto` picks the number of trees: fit at `ntree_tuning.start` trees, multiply by `factor`, and refit with a new seed until two successive fits agree on per-SNP R² (Pearson r), predictor importance (largest change as a fraction of total importance) and turnover space (correlation of pairwise distances between random grid cells in transformed-predictor space). The larger fit of the first agreeing pair is kept. If `ntree_tuning.max` is reached first, that value is used and a warning goes to the `gf_fit` log and `qc.html`. The chosen ntree (`gf/ntree_selected.tsv`) is used for all randomized fits. Set `ntree` to an integer to skip tuning.
* `maxLevel: auto` uses `floor(log2(0.368 * nSites / 2))`, which is 1 for 14 sites.
* The randomizations shuffle whole rows of the predictor table, which keeps the correlations among predictors. Observed and randomized fits use the same `ntree` and `maxLevel`. Empirical p-values are `(1 + #null ≥ obs) / (1 + N)`.
* Turnover colours follow the original script: an unscaled PCA of the transformed predictors, with RGB = (PC1+PC2, −PC2, PC3+PC2−PC1) rescaled to 0–255. Each PC's sign is fixed so that its largest loading is positive.
