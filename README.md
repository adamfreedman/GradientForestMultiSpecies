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
| `popmap` | header `sampleid popid`, columns separated by spaces and/or tabs (or commas): VCF sample → population. Empty fields are an error |
| `latlongs` | CSV, header `popid,sitename,lat,long` |
| `id_regex` | optional column: regex deleted from VCF sample IDs before matching `popmap` (e.g. `^MCZA_`); omit it when the IDs already match |
| `min_ind_per_pop` | optional column: minimum individuals for a population to be kept for this species; blank or absent = `genotypes.min_ind_per_pop` in `config.yaml` |

`config/config.yaml` sets the environmental-layer directory, the study area
(any vector file, or a bounding box), predictor pruning, minimum individuals per population and GF
parameters. Each option is commented in the file.

## Steps

| rule | language | output (under `results/`) |
|---|---|---|
| `study_area` | Python | `env/study_area.geojson`, `env/study_area_bbox.tsv` |
| `get_boundaries` | Python | `env/boundaries.geojson`: country borders for the map (Natural Earth 1:10m by default, downloaded once) |
| `prep_env_stack` | Python | `env/env_stack.tif` (clipped multiband stack, bands named `bio1..bio19`) |
| `validate_inputs` | Python | `{sp}/qc/validation.txt`, `{sp}/qc/sample_pop.tsv`, `{sp}/qc/unmapped_samples.tsv`. Individuals that do not map to a population with coordinates (in the VCF but not the popmap, in the popmap but not the VCF, or in a popid missing from latlongs) are excluded with a warning and listed in `unmapped_samples.tsv` and `qc.html`. Fails, naming the population(s) in the log, if a latlongs population has no individuals in the VCF; also fails on duplicate IDs or sites off the environmental grid |
| `extract_env` | Python | `{sp}/env_sites.tsv` (environment at each population) |
| `vcf_to_popfreq` | Python | `{sp}/genotypes/popfreq.tsv` (populations × SNPs, minor-allele frequency); no MAF or missingness filter; drops only SNPs with no called genotype in some population or no frequency variation among populations (counts in `snp_filter_summary.tsv`); `snp_info.tsv` maps each SNP id (`snp_<VCF ID>`, or `<CHROM>_<POS>` when the VCF has no IDs) to its `CHROM`, `POS` and `VCF_ID` |
| `genotype_pca`, `plot_genotype_pca` | Python, R | `{sp}/pca/`: PCA of individual genotypes, for population structure and sample QC (not used by gradient forest). `variance_explained.tsv` (per PC: eigenvalue, % of variance, cumulative %); `sample_scores.tsv` (one row per individual: `sample_id`, `popid`, `sitename`, `missing_frac`, `in_gf_model`, then its coordinate on every PC); `pca_info.tsv` (individuals and SNPs used); `population_check.tsv` (per individual: distance to its own population's median and to the nearest other population's, and whether it is flagged as sitting with another population; settings in `pca.population_check` in `config.yaml`); `pca_pc1_pc2.pdf/png` and `pca_pc1_pc3.pdf/png` (each population a distinct colour and symbol, with a legend; population IDs at the population medians; flagged individuals ringed in black). Uses every individual mapped to a population with coordinates, including populations too small for the GF model, and every biallelic SNP that varies among them (no MAF or missingness filter); SNPs are scaled by sqrt(2p(1-p)) and missing genotypes set to the mean. Shown in `qc.html`, with a table of flagged individuals |
| `admixture_input`, `admixture_run`, `admixture_pong`, `admixture_summary`, `plot_admixture` | plink2, ADMIXTURE, pong, Python, R | `{sp}/admixture/`: ADMIXTURE ancestry proportions (population structure; not used by gradient forest), for the PCA's individuals and SNPs (no MAF, missingness or LD filter). ADMIXTURE runs `admixture.replicates` times (seeds `seed`, `seed + 1`, ...) at every K from `k_min` to `k_max`, with `cv_folds`-fold cross-validation (`runs/K{k}/rep{r}/`: `.Q`, `.P`, log). `cv_error.tsv` (per run: CV error, log-likelihood); `k_summary.tsv` (per K: mean, SD and minimum CV error, number of pong modes, runs in the largest mode, `stable`); `selected_k.tsv` (the K with the lowest mean CV error, whether it is stable, and the same for K − 1); `ancestry_K{k}.tsv` (per K: one row per individual, `sample_id`, `popid`, then `cluster1..K`, from the representative run of pong's largest mode, with clusters matched across K); `selected_k_modes.tsv` (at the selected K, the representative run of every pong mode); `pong/` (pong's inputs `filemap.txt`, `ind2pop.txt`, `pop_order.txt`, and its results in `pong/output/`, run without the web server or verbose output; to browse the runs interactively, run `pong -m filemap.txt -i ind2pop.txt -n pop_order.txt -o <new dir>` in `pong/`); plots `admixture_selected_k` (selected K; when it is not stable, K − 1 below it, and when K − 1 is not stable either, also the largest stable K below the selected K), `admixture_all_k`, `selected_k_modes` (the 6 largest modes) and `cv_error` (`.pdf/.png`). K is stable when pong's largest mode holds at least `stable_min_frac` of the runs (default 0.9). Individuals are grouped by population (labels as in the PCA plots), with populations of similar ancestry next to each other (average-linkage clustering of the population mean proportions at every K from 2 to the selected K). Shown in `qc.html`, with a warning when the selected K is not stable |
| `install_gf` | R | one-time install of `extendedForest`/`gradientForest` from R-Forge |
| `select_predictors` | R | `{sp}/predictors/candidates.txt` (layers that vary across sites), `selected.txt`, `report.tsv` (rank, screening importance, reason for each drop), `correlation.png`. Of variables correlated above `predictors.corr_threshold`, keeps the one ranked higher by `predictors.rank_by` (default: R²-weighted importance in a screening fit with all candidates) |
| `gf_fit` | R | `{sp}/gf/gf_model.rds`, `importance.tsv`, `snp_rsq.tsv`, `summary.tsv` |
| `gf_randomize` | R | `{sp}/gf/random/batch{n}.tsv`: `n_random` fits (default 1000) with environment rows permuted and the same GF settings, run `random_batch_size` per SLURM job. `batch{n}_snp_rsq.tsv.gz` keeps each replicate's per-SNP R² for later per-SNP p-values |
| `gf_random_summary` | R | `{sp}/gf/random_null.tsv`, `random_pvalues.tsv` |
| `gf_snp_pvalues` | R | `{sp}/gf/snp_pvalues.tsv`: per-SNP `snp` id, VCF `CHROM`, `POS` and `VCF_ID`, observed R², empirical p = (1 + no. randomized fits with R² ≥ observed) / (n_random + 1) and `p_bh`, the Benjamini-Hochberg adjustment (`p.adjust(method = "BH")`). Only SNPs in the observed model (R² > 0) are tested, and BH is over those; SNPs with R² ≤ 0 are listed with `NA` p-values |
| `screen_leave_one_out`, `predictor_stability` | R | `{sp}/predictors/stability_variables.tsv` (per variable: how often it is selected when one population is left out, rank range, what replaces it), `stability_sites.tsv` (per left-out population: variables added or lost). Diagnostic only, shown in `qc.html`; `predictors.stability: false` skips it |
| `gf_predict` | R | `{sp}/gf/turnover_stack.tif` (transformed predictors), `turnover_rgb.tif`, PCA objects |
| `qc_report` | R | `{sp}/qc.html` (`{sp}/refit_sigvars/qc_refit.html` for the refit), one self-contained file with figures embedded: warnings (including ntree not converging), ntree selection, model and randomization summaries, figures, predictor and SNP filtering, genotype PCA and ADMIXTURE (`qc.html` only), populations |
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
* Predictor selection: candidates are ranked (by default by R²-weighted importance in a screening gradient forest with all candidates, `screen_ntree` trees), then taken in rank order, skipping any whose |r| at the sites with an already-kept variable exceeds `corr_threshold`. So of a correlated set the variable that fits the allele frequencies better is kept, and no two kept variables are correlated above the threshold. `keep_priority` variables are ranked first; `rank_by: correlation` ranks by lowest mean |r| instead and ignores the SNPs.
* Sensitivity to sampling design: the selection is repeated once per population with that population left out (screening fit and correlations recomputed). `qc.html` reports how often each selected variable survives and what stands in for it when it does not, so a variable that is often swapped for a correlated one can be read as representing that group. The model itself always uses the full-data selection.
* No selection: `sbatch run_slurm.sh --use-all-predictors` (or `--config use_all_predictors=true`) puts every candidate variable in the model. Results go to `<outdir>_all_predictors`, so a run with selection is never overwritten; the environmental layers are shared.
* Importance floor: `sbatch run_slurm.sh --importance-floor 0.01` (or `--config importance_floor=0.01`) also drops selected variables whose share of the summed screening importance is below 1%, i.e. variables kept only because nothing correlates with them. Default: no floor. Results go to `<outdir>_importance_floor0.01`. Variables in `keep_priority` are exempt; the leave-one-out diagnostic applies the same floor.
* Alternative runs are separate from the default run. With `--use-all-predictors` or `--importance-floor`, the run builds only its own directory (`<outdir>_all_predictors` or `<outdir>_importance_floor<value>`, e.g. `results_importance_floor0.02` for `--importance-floor 0.02`); it does not produce or need the default results, and it leaves existing default results untouched. To get both, launch two runs one after the other (e.g. `sbatch --dependency=afterany:<jobid> run_slurm.sh --importance-floor 0.01`), not at the same time: two Snakemake instances in the same directory interfere with each other. All runs share only the one-time setup under `<outdir>`: the prepared environmental layers and boundary files (`<outdir>/env/`), the gradientForest install marker (`<outdir>/.gf_installed`) and their logs (`<outdir>/logs/`), each built only if missing. Everything else, including the per-species genotype tables, predictor selection, model fits, randomizations and QC reports, is built again in the alternative run's directory.
* Why fit-based selection is the default: [docs/variable_selection.md](docs/variable_selection.md) compares four approaches: (A) the original correlation-only pruning, (B) fit-based selection (the default), (C) all predictors and (D) fit-based selection with a 1% importance floor. It also covers the refit on significant variables, and reports the seed and leave-one-population-out stability tests on the chimpanzee and Trachylepis data.
* The randomizations shuffle whole rows of the predictor table, which keeps the correlations among predictors. Each randomized fit uses the same predictors as the observed model, so the tests ask whether shuffled environments could give the observed fit (per-SNP R², whole-model statistics, importance) under that model. Because the predictors were chosen with the real data, the whole-model and importance p-values are conditional on that choice and somewhat optimistic; the per-SNP tests are conditional on the model by design. Observed and randomized fits use the same `ntree` and `maxLevel`. Empirical p-values are `(1 + #null ≥ obs) / (1 + N)`.
* Turnover colours follow the original script: an unscaled PCA of the transformed predictors, with RGB = (PC1+PC2, −PC2, PC3+PC2−PC1) rescaled to 0–255. Each PC's sign is fixed so that its largest loading is positive.

## Tests

GitHub Actions (`.github/workflows/tests.yml`) runs two jobs on every push and pull request:

* **Unit tests and dry runs** (a few minutes): `pytest` tests of the Python scripts (`tests/test_*.py`) and Snakemake dry runs with the test configuration, which check that the DAG builds, that alternative runs write to their own directories, and that bad settings are rejected. `snakemake --lint` also runs, as advice only.
* **Full workflow on test data** (about 15–25 minutes; less when the cached conda environments can be reused): the whole workflow on a small synthetic dataset, three times (default selection, `importance_floor=0.07`, `use_all_predictors=true`), followed by checks on the outputs (`tests/test_outputs.py`) and R tests of the selection code (`tests/test_gf_common.R`). The QC reports and logs are kept as a downloadable artifact of each run.

The test data in `.test/data` are two simulated species on a 4 × 3 degree grid with six environmental layers, built to exercise the edge cases: correlated, constant and noise layers, a site on a NA cell that has to be snapped, a population below `min_ind_per_pop`, VCF IDs that need `id_regex`, an unmapped sample, a multiallelic site, plain and gzipped VCFs. `.test/make_test_data.py` regenerates them. `.test/config/config.yaml` uses small settings and writes to `.test/results`, so test runs never touch real results.

To run the tests yourself, from the repository root:

```bash
mamba env create -f tests/environment.yaml && mamba activate gfms-test
pytest tests --ignore=tests/test_outputs.py
snakemake --use-conda --cores 4 --configfile .test/config/config.yaml
pytest tests/test_outputs.py
$(ls .snakemake/conda/*/bin/Rscript | head -1) tests/test_gf_common.R
```
