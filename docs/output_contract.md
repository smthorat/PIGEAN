# PIGEAN Engine Output Contract

## Purpose

This document specifies the exact schema of every column produced by the PIGEAN
engine (`priors.py`) across its four output files. Each column is classified by
verification status and display policy, which determines whether it may be shown
to scientists with an interpretation or must be treated as a raw/internal value.

The wrapper (`pigean/interpretation.py`, `pigean/report.py`) uses this contract
to decide what to display and how to describe it. **No metric may be displayed
as an interpreted result unless it is marked VERIFIED + INTERPRET below.**

---

## Classification System

### Verification Status

| Status | Meaning |
|---|---|
| **VERIFIED** | Column source has been traced through `priors.py`, its mathematical definition confirmed, and its scale/units documented. Safe for scientist-facing interpretation with an accompanying interpretation contract. |
| **UNVERIFIED** | Column exists but its interpretation is nuanced, internal-only, or not fully traced. Must not be presented as an interpreted result. |

### Display Policy

| Policy | Meaning |
|---|---|
| **INTERPRET** | Display in the report with an interpretation contract (definition, scale, direction, anti-interpretation). Only for VERIFIED columns. |
| **RAW_ONLY** | May be displayed as a raw value in advanced/debug sections, but without interpretation text. |
| **OMIT** | Do not display in any user-facing report section. Internal, redundant, or confusing for end users. |

### Stochasticity

Columns marked **Stochastic: Yes** will produce different values across runs of
the same input because the engine's Gibbs sampler uses NumPy's random number
generator, which is **not seeded** (`np.random.seed()` is never called in
`priors.py`). Python's `random.seed(0)` is called but this does not affect NumPy.

---

## 1. gs.out — Gene Statistics

One row per gene (~42,216 rows for standard reference). Rows are written sorted
by `combined_prior_Ys` descending (or by `priors` / `Y` as fallbacks). The set
of columns varies by analysis mode.

### Core Columns (present in all modes)

| Column | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|
| `Gene` | Gene identifier from reference gene map | string | No | VERIFIED | INTERPRET |
| `prior` | Gene-level prior log-odds from gene-set model: `X_orig · (betas / scale_factors)`, averaged across Gibbs chains and post-burn-in iterations, mean-centered | log-odds | Yes | VERIFIED | INTERPRET |
| `prior_adj` | `prior` regressed on gene N (number of gene-set memberships) to remove confounding by gene-set size: `prior − slope·N − intercept` | log-odds (adjusted) | Yes | UNVERIFIED | RAW_ONLY |
| `combined` | `prior + log_bf`: combined evidence score integrating gene-set annotations and observed data | log-odds | Yes | VERIFIED | INTERPRET |
| `combined_adj` | `combined` regressed on gene N, analogous to `prior_adj` | log-odds (adjusted) | Yes | UNVERIFIED | RAW_ONLY |
| `combined_D` | Model-derived posterior probability of gene-disease association: `exp(combined) / (1 + exp(combined))` | [0, 1] probability | Yes | VERIFIED | INTERPRET |
| `log_bf` | Gene-level log Bayes factor from observed data. In positive-controls mode this is the input signal itself; in GWAS/exome modes it is engine-derived | log-odds | Mode-dependent | VERIFIED | INTERPRET |
| `N` | Number of gene sets containing this gene (sum of absolute `X_orig` weights plus ignored gene-set contributions) | count | No | VERIFIED | INTERPRET |
| `Chrom` | Chromosome from gene location reference | string | No | VERIFIED | INTERPRET |
| `Start` | Genomic start position (base pairs) | integer | No | VERIFIED | INTERPRET |
| `End` | Genomic end position (base pairs) | integer | No | VERIFIED | INTERPRET |

### Mode-Specific Columns

| Column | Present in | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|---|
| `positive_control` | positive-controls mode | Input echo of user-provided positive-control annotation value | float (score) | No | VERIFIED | INTERPRET |
| `huge_score` | GWAS+exome combined | Sum of `huge_score_gwas` + `huge_score_exomes` | log-odds | No | UNVERIFIED | RAW_ONLY |
| `huge_score_gwas` | GWAS mode | GWAS signal-to-gene (s2g) score: log Bayes factor from GWAS proximity mapping | log-odds | No | UNVERIFIED | RAW_ONLY |
| `huge_score_gwas_uncorrected` | conditional | GWAS huge score before covariate correction | log-odds | No | UNVERIFIED | OMIT |
| `huge_score_exomes` | exome mode | Exome-data signal-to-gene score | log-odds | No | UNVERIFIED | RAW_ONLY |

### Conditional Columns (appear based on engine configuration)

| Column | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|
| `combined_se` | Standard error of combined prior+Y estimate (often None / not computed) | log-odds SE | Yes | UNVERIFIED | OMIT |
| `log_bf_regression` | Log BF used for regression; may differ from `log_bf` if adjusted | log-odds | Yes | UNVERIFIED | OMIT |
| `log_bf_uncorrected` | Y before covariate correction via `_correct_huge` | log-odds | No | UNVERIFIED | OMIT |
| `log_bf_w` | Whitened log BF: `L⁻¹ Y` where L is the Cholesky factor of gene correlation matrix | whitened log-odds | No | UNVERIFIED | OMIT |
| `log_bf_fw` | Fully whitened: `(LL')⁻¹ Y` for GLS-style inference | inverse-correlation-weighted | No | UNVERIFIED | OMIT |
| `prior_orig` | Prior from before the most recent Gibbs update (snapshot) | log-odds | Yes | UNVERIFIED | OMIT |
| `prior_adj_orig` | Adjusted prior from before the most recent update | log-odds | Yes | UNVERIFIED | OMIT |
| `batch` | Batch assignment for the gene | string | No | UNVERIFIED | OMIT |
| *`<covariate_name>`* | Z-scored gene-level covariates (dynamic, one column per covariate excluding intercept) | Z-score | No | UNVERIFIED | OMIT |

---

## 2. gss.out — Gene-Set Statistics

One row per gene set (~22,000–42,000 rows depending on mode and filtering). Rows
sorted by `betas / scale_factors` descending.

### Standard Columns (26 observed in typical runs)

| Column | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|
| `Gene_Set` | Gene-set identifier | string | No | VERIFIED | INTERPRET |
| `label` | Source annotation file label | string | No | VERIFIED | INTERPRET |
| `N` | Gene-set size: total absolute weight of genes in this set (sum of `|X_orig|` column) | count | No | VERIFIED | INTERPRET |
| `scale` | Scale factor used to normalize X columns; dividing internal beta by this gives external (per-gene) beta | unitless | No | UNVERIFIED | RAW_ONLY |
| `beta_tilde` | Marginal OLS/logistic regression effect estimate in external (per-gene) units | per-gene effect | Yes | UNVERIFIED | RAW_ONLY |
| `beta_tilde_internal` | Same in internal (scaled) units | internal | Yes | UNVERIFIED | OMIT |
| `P` | Marginal regression p-value: `2 · Φ(−|Z|)` where `Z = beta_tilde / SE` | [0, 1] | Yes | VERIFIED | INTERPRET |
| `Z` | Z-score from marginal regression: `beta_tilde / SE` | Z-score | Yes | UNVERIFIED | RAW_ONLY |
| `SE` | Standard error of `beta_tilde` in external units | per-gene SE | Yes | UNVERIFIED | RAW_ONLY |
| `beta` | Posterior mean effect from spike-and-slab sampler, corrected for LD between gene sets (V matrix), in external units. **Primary gene-set effect estimate.** | per-gene log-odds contribution | Yes | VERIFIED | INTERPRET |
| `beta_internal` | Same in internal (scaled) units | internal | Yes | UNVERIFIED | OMIT |
| `beta_uncorrected` | Posterior mean assuming independent gene sets (no V correction), external units | per-gene effect | Yes | UNVERIFIED | RAW_ONLY |
| `avg_postp` | Posterior inclusion probability (PIP) for this gene set: probability of non-zero effect | [0, 1] | Yes | VERIFIED | INTERPRET |
| `beta_tilde_orig` | Original `beta_tilde` before Gibbs re-estimation | per-gene effect | Yes | UNVERIFIED | OMIT |
| `beta_tilde_internal_orig` | Same in internal units | internal | Yes | UNVERIFIED | OMIT |
| `P_orig` | Original p-value before Gibbs | [0, 1] | Yes | UNVERIFIED | OMIT |
| `Z_orig` | Original Z-score | Z-score | Yes | UNVERIFIED | OMIT |
| `SE_orig` | Original SE | per-gene SE | Yes | UNVERIFIED | OMIT |
| `beta_orig` | Original posterior mean beta | per-gene effect | Yes | UNVERIFIED | OMIT |
| `beta_internal_orig` | Same in internal units | internal | Yes | UNVERIFIED | OMIT |
| `beta_uncorrected_orig` | Original uncorrected beta, external | per-gene effect | Yes | UNVERIFIED | OMIT |
| `beta_uncorrected_internal_orig` | Same in internal units | internal | Yes | UNVERIFIED | OMIT |
| `avg_cond_beta_orig` | Original average beta conditional on non-zero: `beta / postp` | per-gene effect | Yes | UNVERIFIED | OMIT |
| `avg_postp_orig` | Original PIP | [0, 1] | Yes | UNVERIFIED | OMIT |
| `p_used` | Prior probability of inclusion (spike-and-slab p) used for this gene set | [0, 1] | Yes | UNVERIFIED | RAW_ONLY |
| `sigma2_used` | Prior variance of beta for this gene set: `sigma2 · scale_factor^sigma_power` (no thresholding) | variance | Yes | UNVERIFIED | RAW_ONLY |

### Additional Conditional Columns

| Column | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|
| `inf_beta` | Infinitesimal (ridge-like) posterior mean beta, external units | per-gene effect | No | UNVERIFIED | RAW_ONLY |
| `inf_beta_orig` | Original infinitesimal beta | per-gene effect | No | UNVERIFIED | OMIT |
| `avg_cond_beta` | Average beta conditional on non-zero: `beta / postp` (from `_calculate_non_inf_betas`, not Gibbs) | per-gene effect | Yes | UNVERIFIED | RAW_ONLY |
| `sigma2_thresholded` | Prior variance with sigmoid soft-thresholding applied | variance | Yes | UNVERIFIED | OMIT |
| `O` | Overlap score: sum of squared pairwise V correlations (measures collinearity with other gene sets) | sum of squared correlations | No | UNVERIFIED | OMIT |
| `X_O` | Transformed overlap: `O / SE²`, further scaled by `scale_factors^sigma_power` | scaled overlap | No | UNVERIFIED | OMIT |
| `weight` | Overlap-derived weight used in sampling | weight | No | UNVERIFIED | OMIT |
| *`avg_<covariate>`* | Average value of each gene-level QC covariate across genes in this set (Z-scored), one column per covariate excluding intercept, plus `avg_huge_adjustment` | Z-score or adjustment | No | UNVERIFIED | OMIT |
| `avg_avg_metric` | Mean of the mean QC metrics across genes in this set | Z-score (mean of means) | No | UNVERIFIED | OMIT |

---

## 3. ggss.out — Gene-Gene-Set Statistics

Sparse long-format file. Each row is one (gene, gene_set) pair where the gene
belongs to the gene set AND the gene set's beta passes the output filter. Rows
sorted by gene (descending `combined_prior_Ys`) then by gene set within gene
(descending `beta / scale_factor`). Typical size: ~200–284,000 rows.

### Core Columns (7)

| Column | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|
| `Gene` | Gene identifier | string | No | VERIFIED | INTERPRET |
| `prior` | Gene-level prior log-odds (same value as in gs.out) | log-odds | Yes | VERIFIED | INTERPRET |
| `combined` | Combined prior + data (same value as in gs.out) | log-odds | Yes | VERIFIED | INTERPRET |
| `log_bf` | Gene log Bayes factor (same value as in gs.out) | log-odds | Mode-dependent | VERIFIED | INTERPRET |
| `gene_set` | Gene-set identifier | string | No | VERIFIED | INTERPRET |
| `beta` | Gene-set beta in external units (prefers corrected betas; falls back to `beta_tildes`) | per-gene effect | Yes | VERIFIED | INTERPRET |
| `weight` | Weight of this gene in this gene set from the X matrix (membership weight, often 0 or 1) | raw weight | No | VERIFIED | INTERPRET |

### Conditional Columns

| Column | Present when | Source | Scale | Stochastic | Verified | Display |
|---|---|---|---|---|---|---|
| `log_bf_for_regression` | when `Y ≠ Y_for_regression` | Log BF used for regression (may be adjusted) | log-odds | Yes | UNVERIFIED | OMIT |
| `huge_score_gwas` | when both GWAS + exomes exist | GWAS huge score (dynamically chosen label) | log-odds | No | UNVERIFIED | RAW_ONLY |
| `huge_score` | when only one of GWAS/exomes exists | Whichever huge score exists | log-odds | No | UNVERIFIED | RAW_ONLY |
| `huge_score_exomes` | when both GWAS + exomes exist | Exome huge score | log-odds | No | UNVERIFIED | RAW_ONLY |

---

## 4. p.out — Model Parameters

Three-column TSV: `Parameter`, `Version`, `Value`. Parameters that change during
the run (e.g., due to Gibbs restarts or hyperparameter updates) appear as
multiple rows with incrementing `Version` numbers (1-indexed). Most parameters
have a single version.

### Convergence-Relevant Parameters

| Parameter | Type | Meaning | Notes |
|---|---|---|---|
| `num_gibbs_iter` | int | Final Gibbs iteration number, **0-indexed**. Value 499 means 500 iterations were executed (loop `range(500)` produces 0…499). If this equals `max_num_iter − 1`, the iteration cap was hit. | **Primary convergence indicator.** |
| `num_gibbs_restarts` | int | Number of Gibbs restarts. Value > 0 indicates the sampler had difficulty converging and hyperparameters were adjusted. | **Primary convergence indicator.** |
| `num_chains` | int | Number of outer Gibbs MCMC chains (default 10). | Informational. |
| `p` | float, multi-version | Prior inclusion probability (spike-and-slab). Multiple versions arise from hyperparameter adjustments during restarts. | Trajectory of model fitting. |
| `sigma2` | float, multi-version | Prior variance of gene-set effects. Multiple versions from restarts. | Trajectory of model fitting. |

### Informational Parameters

| Parameter | Type | Meaning |
|---|---|---|
| `num_chains_betas` | int | Number of chains for inner beta sampling (default 5) |
| `sigma2_cond` | float | `sigma2 / p` (conditional prior variance) |
| `sigma_power` | int | Power for scale-dependent sigma |
| `num_gene_sets_read` | int | Final count of gene sets loaded after filtering |
| `num_genes_read` | int | Final count of genes loaded |
| `mad_threshold` | int | MAD threshold for outlier chain detection (default 10) |
| `use_mean_betas` | bool | Whether mean betas were used |
| `sparse_solution` | bool | Whether sparse solution was used |
| `sparse_frac` | float | Fraction for sparse solution |
| `sparse_max` | float | Maximum for sparse solution |
| `sparse_frac_betas` | float | Fraction for sparse betas |
| `filter_gene_set_p` | float | P-value threshold for gene-set filtering at read time |
| `filter_negative` | bool | Whether negative-weight gene sets were filtered |
| `threshold_weights` | float | Weight threshold for binarizing X |
| `cap_weights` | bool | Whether gene-set weights were capped |
| `max_num_gene_sets` | int | Maximum number of gene sets |
| `filter_gene_set_metric_z` | float | Z-score threshold for gene-set metric filtering |
| `sigma_num_devs_to_top` | float | Number of std devs for sigma top estimation |
| `p_noninf_inflate` | float | Inflation factor for `p_noninf` |
| `num_X_batches` | int | Number of batches for X processing |
| `read_X_run_logistic` | bool | Whether logistic regression was used in X processing |
| `gene_set_prune_threshold` | float | Threshold for pruning correlated gene sets |
| `gene_set_prune_deterinistically` | str/None | Whether pruning is deterministic |
| `max_allowed_batch_correlation` | float | Correlation threshold for batch assignment |
| `initial_linear_filter` | bool | Whether initial linear filter was applied |
| `correct_betas_mean` | bool | Whether beta mean correction was applied |
| `correct_betas_var` | bool | Whether beta variance correction was applied |
| `adjust_priors` | bool | Whether priors were adjusted post-Gibbs |

### Restart-Related Parameters (appear only when restarts > 0)

| Parameter | Type | Meaning |
|---|---|---|
| `p_adj` | float, multi-version | Adjusted p value after restart |
| `sigma2_adj` | float, multi-version | Adjusted sigma2 after restart |
| `p_scale_factor` | float | Scale factor applied to p during restart |
| `fraction_required_to_not_increase_hyper` | float | Fraction threshold that triggered restart |

---

## 5. Convergence Assessment from p.out

The wrapper assesses convergence using two parameters from p.out:

1. **`num_gibbs_iter`**: If this equals `max_num_iter − 1` (default: 499 for
   max 500), the sampler hit its iteration cap. This is a convergence concern
   but does **not** definitively mean non-convergence — the SEM criterion may
   still have been met (confirmed via log file parsing).

2. **`num_gibbs_restarts`**: If > 0, the sampler required hyperparameter
   adjustments to make progress. This indicates initial convergence difficulty.

The log file provides additional evidence:
- Final SEM ratio (`max_ratio` in the log) compared to the 0.01 threshold
- R-hat convergence statistics
- The rare "Desired Gibbs precision achieved" message (requires non-default
  `--increase-hyper-if-betas-below` option to be set)

See `pigean/convergence.py` for the full two-tier assessment logic.

---

## 6. Stochasticity Note

All columns marked **Stochastic: Yes** will produce different numerical values
across runs of the same input data. This is because the engine's Gibbs sampler
relies on `numpy.random`, which is **not explicitly seeded** in `priors.py`.
While `random.seed(0)` is called for Python's built-in `random` module, this
does not affect NumPy's random state.

Consequently:
- Gene rankings by `combined_D` may shift between runs
- Gene-set `beta` and `avg_postp` values will differ
- Row counts in `ggss.out` may vary (stochastic filtering of gene sets)
- Parameter versions in `p.out` for `p` and `sigma2` will differ

Deterministic columns (identifiers, coordinates, input echoes, counts) remain
stable across runs.

---

## 7. Column Presence by Analysis Mode

| Column | positive-controls | bf | z-score | percentile | exome | gwas |
|---|---|---|---|---|---|---|
| `positive_control` (gs) | ✓ | — | — | — | — | — |
| `huge_score_exomes` (gs) | — | — | — | — | ✓ | — |
| `huge_score_gwas` (gs) | — | — | — | — | — | ✓ |
| `huge_score` (gs) | — | — | — | — | ✓+gwas | ✓+exome |

All other core columns are present in all modes. The gss.out and p.out schemas
are identical across modes. The ggss.out schema is identical except for the
optional `huge_score_*` columns.
