# Advanced PIGEAN Modes

## Two independent choices

`--analysis` selects the evidence adapter (the data a researcher supplies).
`--mode` selects the statistical pathway used by the frozen engine. The wrapper
resolves both before execution and never falls back to another mode.

## Compatibility matrix

| Evidence type | `standard` | `naive-priors` |
|---|---|---|
| positive-controls | SUPPORTED | SUPPORTED |
| gene-bayes-factor | SUPPORTED | NOT_VERIFIED |
| exome | SUPPORTED | NOT_VERIFIED |
| gwas | SUPPORTED | NOT_VERIFIED |
| gene-z-score | ENGINE_BLOCKED | ENGINE_BLOCKED |
| gene-percentile | ENGINE_BLOCKED | ENGINE_BLOCKED |

All four gene-set profiles (`default`, `mouse-only`, `msigdb-only`, and
`custom`) remain available to both exposed modes. Background-list handling uses
the same positive-control adapter path. Only the default positive-control
profile was exercised end-to-end for naive-priors in Phase 6; alternative
profiles/backgrounds inherit validated preprocessing but are not separate
scientific benchmarks.

## `standard`

- Engine subcommand: `gibbs`
- Engine function: `GeneSetData.run_gibbs`
- Meaning: the Phase 1–5 validated outer gene-prior Gibbs workflow
- Required inputs: those of the selected supported evidence adapter, gene-set
  matrices, gene map, and genome-build references
- Advanced parameters: none exposed
- Outputs: standard `gs.out`, `gss.out`, `ggss.out`, and `p.out`
- Interpretation: FULLY_COMPATIBLE with the Phase 5 contract
- PIGEAN stability: APPLICABLE
- Formal MCMC convergence: NOT_FORMALLY_ASSESSED
- Status: SUPPORTED

## `naive-priors`

- Engine subcommand: `naive_priors` (singular `naive_prior` is an engine alias,
  but the wrapper emits only the canonical token)
- Engine functions: `calculate_gene_set_statistics`,
  `calculate_non_inf_betas`, and `calculate_naive_priors`
- Meaning: estimates conditional gene-set effects with the engine's inner
  sampler, computes `X_orig · (betas / scale_factors)`, mean-centers across
  present/missing genes, and adds input `Y` to form `combined`
- Required input: positive controls plus the same gene-set/gene-map/reference
  inputs as standard mode
- Advanced parameters: none exposed
- Outputs: all four core files; `gs.out` omits `combined_D` and adjusted
  columns, and `p.out` has no outer-Gibbs iteration/restart fields
- Interpretation: PARTIALLY_COMPATIBLE. Gene ranking uses `combined`, and the
  prior is explicitly described as the naive estimator.
- PIGEAN stability: NOT_APPLICABLE to the bypassed outer Gibbs loop
- Formal MCMC convergence: NOT_APPLICABLE to the bypassed outer Gibbs loop
- Status: SUPPORTED for positive controls only

Example:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --mode naive-priors \
  --input examples/gene_list \
  --output results/runs/example_naive
```

The mode is faster in the recorded Phase 6 control run, but the estimators and
schemas differ. The comparison is descriptive and does not establish a
scientifically preferred mode.

## Investigated but not exposed

### Factor / `naive_factor` — NOT_VERIFIED

The engine performs Bayesian non-negative matrix factorization with automatic
relevance determination over a gene-set × gene or gene-set × phenotype matrix.
Its input matrix, masks, output loadings, clusters, and relevance scores change
with nine documented anchor scenarios. It can also project factors across
PheWAS matrices. These are specialized workflows, not a single safe extension
of the wrapper's `gs/gss/ggss` interpretation contract. `factor` and
`naive_factor` are therefore not registered.

### PheWAS — NOT_VERIFIED

The engine reads long-form gene × phenotype values (`Gene`, `Pheno`, and one
or more of `log_bf`, `combined`, `prior`) and/or gene-set × phenotype effect
matrices. It runs phenotype-wise regressions and can feed factorization. The
units are genes/gene sets by phenotypes, phenotypes are processed in batches,
and optional factor projection makes output semantics mode-dependent. There is
no Phase 5 adapter or interpretation contract for these outputs, so PheWAS is
not exposed.

### Anchors — INTERNAL_ONLY

`--anchor-phenos`, `--anchor-any-pheno`, `--anchor-genes`,
`--anchor-any-gene`, and `--anchor-gene-set` construct masks/relevance weights
for different factorization matrices. They require specific paired PheWAS
matrices in most cases, and precedence rules can ignore other supplied inputs.
They are not general labels or generic covariates and are withheld.

### `phi` — INTERNAL_ONLY

Engine default: `0.05` (float). In the active Bayesian NMF calculation it is
multiplied by the variance of the non-negative matrix and scales L2/ARD
regularization; larger values are described by the engine as yielding fewer
factors. No scientifically defensible wrapper range or cross-matrix scaling
contract is established.

### `alpha0` — INTERNAL_ONLY

Engine default: `10` (float). It enters the inverse-gamma/ARD relevance prior
through `C = (N + M)/2 + alpha0 + 1` and the internally derived `b0`. A safe
researcher-facing range and sensitivity contract have not been established.

### `beta0` — ENGINE_BROKEN / NOT EXPOSED

Engine default: `1` (float). It is parsed and passed to `run_factor`, but the
active `_bayes_nmf_l2_extension` call does not consume it; `b0` is derived from
`alpha0` instead. Exposing it would falsely imply that it changes the fitted
model.

## Validation behavior

Unknown modes, not-verified evidence combinations, convergence traces for
naive-priors, and internal knobs (`anchor`, `phi`, `alpha0`) fail before engine
execution. The failure manifest records `validation.status = FAIL` and
`execution.status = NOT_RUN`.
