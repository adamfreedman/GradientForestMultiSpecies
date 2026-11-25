# How the predictor-selection method was chosen

This document records the analyses behind the workflow's choice of environmental
predictors for gradient forest (GF). It covers what was tried, the results for
the two test species, and why the current default was adopted. It is written
so the reasoning can be lifted into a methods section or supplement.

Test data: chimpanzee (16 populations) and *Trachylepis affinis* (14
populations), each with the 19 WorldClim bioclimatic variables (bio1–bio19)
extracted at the population sites. All models used 500 trees, which the ntree
tuning step chose for both species. Significance came from 1,000 randomizations
that permute the environmental rows across populations. Per-SNP tests and
Benjamini–Hochberg (BH) correction cover only SNPs with R² > 0.

---

## 1. The problem

The 19 bioclimatic variables are strongly correlated across a few sites. Many
pairs have |r| > 0.9 (for example, chimp bio16 and bio19 have r = 0.98). GF
splits importance among correlated predictors, which makes variable-level
results hard to read. Correlated predictors also add noise to the randomization
null. The standard fix is to keep one variable from each group of correlated
variables (|r| > 0.8 at the sites). The open question is **which member of each
group to keep**.

## 2. Methods compared

| Method | How a variable is chosen within each correlated group | Status |
|---|---|---|
| A. Correlation pruning (original) | In each pair with \|r\| > 0.8, drop the variable with the higher mean \|r\| to the remaining candidates; repeat until no pair exceeds 0.8. Uses only the environment, not the genetic data. | Replaced |
| B. Fit-based selection | Fit a screening GF with every candidate variable. Rank variables by R²-weighted importance. Go down the ranking and keep a variable unless it has \|r\| > 0.8 with one already kept. Uses the genetic data. | **Current default** |
| C. All predictors | No selection; all 19 variables. | Optional (`--use-all-predictors`) |
| Refit on significant variables | After any of A–C, refit using only variables whose importance or accuracy beats the randomization null (p ≤ 0.05). | Produced alongside the main model as a secondary analysis |
| Leave-one-population-out (LOO) selection | Repeat B with each population left out in turn. Considered as a possible replacement for B; adopted as a diagnostic. | **Diagnostic in QC report** |

**Why B over A, in one sentence:** method A picks among near-equivalent
variables on a criterion unrelated to the biology, so the choice is essentially
arbitrary. Method B keeps the member of each group that best explains the allele
frequencies, so selection follows the data. The rest of this document tests
whether B's choice is reliable and whether it matters for the results.

## 3. Results

### 3.1 Variables selected

| | Chimpanzee | Trachylepis |
|---|---|---|
| A. Correlation pruning | bio2 bio3 bio4 bio6 bio10 bio15 bio18 bio19 (8) | bio3 bio4 bio5 bio8 bio9 bio14 bio15 bio18 bio19 (9) |
| B. Fit-based | bio2 bio3 bio4 bio5 bio6 bio8 bio15 bio18 bio19 (9) | bio3 bio4 bio5 bio8 bio9 bio12 bio14 bio15 bio18 bio19 (10) |
| Refit after A | bio19 bio4 | bio8 bio18 bio19 bio3 |
| Refit after B | bio19 bio4 | bio12 bio8 bio18 bio19 bio3 |
| Refit after C | bio19 bio4 bio16 | bio18 bio12 bio6 bio8 bio2 bio7 bio13 bio3 bio19 bio16 |

The two methods differ by only one or two variables:

- **Chimp:** B swapped bio10 for bio5 and bio8. All three are low-importance
  temperature variables.
- **Trachylepis:** B added bio12 (annual precipitation), which ranks 3rd of 19
  in the screening fit. Method A had dropped it because it is correlated with
  many other variables.

### 3.2 Model fit and SNP counts

| Model | SNPs with R² > 0 (null mean; p) | Mean R² (p) | SNPs at BH ≤ 0.05 | SNPs at BH ≤ 0.1 |
|---|---|---|---|---|
| **Chimpanzee** | | | | |
| A. Correlation pruning | 2,772 (826; 0.012) | 0.237 (0.004) | 960 | 2,390 |
| B. Fit-based | 2,675 (805; 0.015) | 0.246 (0.007) | 919 | 2,355 |
| C. All predictors | 2,618 (784; 0.014) | 0.245 (0.012) | 785 | 2,363 |
| Refit after A or B (bio19, bio4) | 2,904 (764; 0.007) | 0.284 (0.005) | 1,401 | 2,714 |
| Refit after C (bio19, bio4, bio16) | 2,868 (730; 0.009) | 0.287 (0.005) | 1,490 | 2,761 |
| **Trachylepis** | | | | |
| A. Correlation pruning | 15,701 (2,382; 0.001) | 0.436 (0.001) | 13,939 | 15,537 |
| B. Fit-based | 16,234 (2,431; 0.001) | 0.438 (0.002) | 14,280 | 15,987 |
| C. All predictors | 15,702 (2,318; 0.001) | 0.460 (0.002) | 13,644 | 15,517 |
| Refit after A | 15,365 (1,971; 0.001) | 0.492 (0.001) | 13,971 | 15,318 |
| Refit after B | 16,071 (2,069; 0.001) | 0.482 (0.001) | 14,293 | 16,017 |
| Refit after C | 15,659 (2,183; 0.001) | 0.529 (0.002) | 14,394 | 15,590 |

Overlap of significant SNPs (BH ≤ 0.1) between B and C: chimp 2,134 shared
(91% of B), trachylepis 15,052 shared (94% of B). Refits: chimp 2,553 shared
(94%), trachylepis 14,760 shared (92%).

**Interpretation**

- Every method gives a significant whole-model fit, and similar SNP counts that
  are mostly the same SNPs. **The SNP-level conclusions do not depend on how
  predictors are selected.**
- **Why B found fewer significant SNPs than A in chimp but more in
  trachylepis.** We wanted to know whether this inconsistency pointed to a
  problem with B. It does not:
  - In chimp, B added two low-importance variables, bio5 and bio8. Weak
    predictors raise the randomization null a little, so slightly fewer SNPs
    pass.
  - In trachylepis, B kept bio12, a strong predictor, so slightly more SNPs pass.
  - The differences are small: −1.5% (chimp) and +2.9% (trachylepis) at
    BH ≤ 0.1.
  - The trachylepis signal is very strong (16,000 SNPs against a null of
    ~2,400), which suggests that much of it follows population structure that
    lines up with the environment.
- **All predictors (C):** fitting all 19 variables doesn't improve the fit. It
  does make the variable-level results harder to read. Correlated partners
  share importance and pass the randomization test together: in chimp, bio16
  (r = 0.98 with bio19) joins bio19; in trachylepis, 10 of 19 variables pass,
  against 5 after selection. **Selection matters for interpreting individual
  variables, not for which SNPs are detected.**

### 3.3 Is the fit-based choice stable? Screening-fit seeds

Method B relies on one screening GF fit. To check whether random-forest noise
could change the choice, we repeated the screening fit with 20 random seeds per
species.

- **Rankings are effectively fixed.** Each variable's importance varied by
  1–2% across seeds (coefficient of variation ≤ 0.02). Neighbouring ranks
  rarely traded places, and never in a way that changed the selection.
- **All 20 seeds gave the same selected set as the pipeline**, for both
  species. Selecting on importance averaged over seeds gave the same set again.
- SNPs with R² > 0 in the screening fit: chimp 2,590–2,635; trachylepis
  15,672–15,734.

Conclusion: **method B is not sensitive to the random seed.** Averaging several
screening fits would cost more and change nothing.

### 3.4 Is the fit-based choice stable? Leave one population out

The selection could still depend on *which populations were sampled*: where
they lie in environmental space determines both the correlations and the
importances. To test this, we repeated method B 16 times for chimp and 14 times
for trachylepis, leaving out one population each time. Correlations were
recomputed among the remaining sites.

| | Chimpanzee | Trachylepis |
|---|---|---|
| Leave-one-out fits that reproduce the full-data set exactly | 7 of 16 | 2 of 14 |
| Variables kept in every fit | bio19, bio4, bio2 | bio19, bio3, bio5 |
| Least stable kept variables | bio8 (kept in 56% of fits) | bio14 (43%), bio15 (50%), bio9 (71%) |
| Most influential populations | Mount Cameroon (3 kept variables lost, 2 others gained), Lobeke, Minta | LOPG (3 lost, 2 gained), KRI (3 lost, 1 gained) |

**Interpretation**

- **The strongest predictors are stable.** The variables that pass the
  randomization test and make up the refits (chimp bio19 and bio4; trachylepis
  bio18, bio19, bio12, bio3, bio8) are kept in 86–100% of leave-one-out fits.
- **Instability is confined to weak predictors in tight correlated groups.**
  When a variable drops out, a highly correlated partner almost always replaces
  it:
  - chimp: bio8 is replaced by bio10, bio6, bio5 or bio9 (|r| 0.81–0.91);
  - trachylepis: bio14 and bio15 are replaced by bio17, part of a dry-season
    precipitation group (bio14–bio17 |r| = 0.97). Which member of that group
    is kept depends on the sampled sites, but the group itself is always
    represented.
- One population can noticeably change the model: leaving out Batanga raised
  chimp SNPs with R² > 0 from ~2,600 to 3,306.
- **Using LOO to choose variables was considered and rejected.** Ranking
  variables by how often they are kept, or by their average importance across
  the leave-one-out fits, reproduced exactly the sets chosen by method B, so
  it adds computation without changing the answer. Requiring a variable to be
  kept in ≥ 75% of fits would instead drop the whole dry-season precipitation
  group in trachylepis, leaving that environmental axis unrepresented. So LOO
  is better used to *report* how sensitive the selection is to sampling than
  to *make* the selection.

### 3.5 Why the refit is not the primary model

The refit on significant variables gives the highest mean R² and the most
significant SNPs (chimp 2,714 against 2,355 at BH ≤ 0.1). However, it selects
variables using the same randomization test it is then evaluated with. That
test permutes environmental rows across populations, and with 14–16
populations its power depends on how evenly the populations cover environmental
space. A few influential populations, as the LOO analysis shows, can drive which
variables pass. We therefore keep the full selected-set model (B) as the primary
result and report the refit as a secondary analysis.

## 4. Decision

1. **Default: fit-based selection (method B)** with |r| threshold 0.8 and
   ranking by R²-weighted importance in a 500-tree screening fit. It uses the
   genetic data to choose among near-equivalent correlated variables, it is not
   affected by seed, and it gives SNP results consistent with other methods.
2. **Diagnostic: leave-one-population-out stability**, run automatically and
   reported in the QC report (`predictors/stability_variables.tsv`,
   `predictors/stability_sites.tsv`). It shows how often each variable is kept,
   what replaces it, and which populations change the selection. The QC report
   warns about kept variables retained in < 75% of fits. It can be turned off
   with `predictors.stability: false`.
3. **Option: `--use-all-predictors`** fits all candidate variables with no
   selection, into a separate output directory (`<outdir>_all_predictors`).
   This is a sensitivity check for selection. Without the switch, selection is
   always done.
4. **Secondary: refit on significant variables**, kept as a separate model and
   report (`refit_sigvars/`), not the primary result.

## 5. Known limitation, not addressed

The correlation rule keeps any variable that is not correlated with a
higher-ranked one, however unimportant. Some kept variables explain almost
nothing: trachylepis bio5 (importance 0.0005, rank 19 of 19) and chimp bio18
(rank 19). Weak predictors slightly raise the randomization null (section 3.2).
A minimum-importance threshold would remove them, but it is not implemented
because no principled cut-off has been settled on.

## 6. Draft narrative for a methods section

> Environmental predictors were the 19 WorldClim bioclimatic variables. Because
> many are highly correlated across sampled populations, we retained one
> variable from each group of correlated predictors (|Pearson r| > 0.8 at the
> sampled sites). Rather than choosing among correlated variables on
> environmental criteria alone, we fit a screening gradient forest with all
> candidate variables (500 trees), ranked variables by R²-weighted importance,
> and retained variables in rank order, skipping any correlated (|r| > 0.8) with
> one already retained. Importance rankings were invariant across 20 random
> seeds of the screening fit, and every seed produced the same predictor set.
> To assess sensitivity to the sampling of populations, we repeated selection
> with each population left out in turn. Predictors that were significant in
> randomization tests were retained in 86–100% of leave-one-out selections;
> instability was confined to weak predictors within tightly correlated groups,
> where a correlated variable (|r| ≥ 0.8) took their place. As sensitivity
> analyses, models fit with all 19 predictors or with predictors chosen by
> correlation alone identified largely the same SNPs (91–94% overlap at FDR
> ≤ 0.1) with similar overall fit. Without selection, however, correlated
> variables shared importance, which complicated interpretation of individual
> predictors.
