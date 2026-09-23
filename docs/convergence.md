# PIGEAN Stability and Convergence Documentation

## Overview

The Phase 5 stability assessment evaluates PIGEAN's engine-specific
max-fractional-SEM criterion.  This is **not** formal MCMC convergence.

## Key Terminology

### PIGEAN Stability Criterion

The PIGEAN engine computes a **max fractional SEM** at each Gibbs iteration:

```
max_fractional_sem = max_over_all_genes(
    SEM_of_cross_chain_running_average_log_posterior_odds
) / range_of_gene_average_log_posterior_odds
```

This measures how well the 10 chains agree on each gene's posterior
estimate.  The threshold is `--max-frac-sem = 0.01` (from the original
`priors.py`, line 498).

**This criterion was defined by the PIGEAN authors.  It was NOT introduced
by the Phase 5 wrapper.**

### Status Values

| Status | Meaning |
|--------|---------|
| `PIGEAN_STABILITY_CRITERION_MET` | The engine's max-fractional-SEM criterion was satisfied |
| `PIGEAN_STABILITY_CRITERION_UNCERTAIN` | Could not verify the criterion from available data |
| `PIGEAN_STABILITY_CRITERION_NOT_MET` | The criterion was not satisfied (restarts + cap hit) |
| `NOT_ASSESSED` | Insufficient data to assess |
| `TRACE_NOT_AVAILABLE` | Trace output not available |
| `ENGINE_FAILED` | The engine did not complete |
| `NOT_APPLICABLE` | The selected mode does not run the outer gene-prior Gibbs loop |

### `PIGEAN_STABILITY_CRITERION_MET` Does NOT Mean

- Formal MCMC convergence established
- R-hat < threshold for all parameters
- Effective sample size is adequate
- Chain traces have been verified for stationarity
- Posterior samples are independent

It **does** mean the engine's built-in criterion is satisfied: the maximum
normalized SEM of cross-chain gene posterior estimates is below the
engine's own threshold.

### Formal MCMC Convergence

```
formal_convergence_status = NOT_FORMALLY_ASSESSED
```

**Reason**: The frozen PIGEAN engine does not retain sufficient independent
chain traces for standard multi-chain diagnostics such as R-hat or
effective sample size.

Specific gaps:
1. Individual chain traces for gene-level log-posteriors are not retained
   (only running sums)
2. Per-chain post-burn-in samples are accumulated into running sums;
   individual samples are discarded
3. Autocorrelation information needed for ESS is not computed
4. Chains share identical initialization (not overdispersed)

### Naive-priors mode

The `naive-priors` engine route bypasses `run_gibbs`, the outer gene-prior
Gibbs loop assessed by the PIGEAN max-fractional-SEM criterion. For this mode:

- PIGEAN stability: `NOT_APPLICABLE`
- Formal MCMC convergence: `NOT_APPLICABLE`
- Iteration status: `NOT_APPLICABLE`

The mode still calls the stochastic inner gene-set effect sampler. That does
not make the outer-Gibbs diagnostic applicable and is stated separately in the
report caveats.

## Assessment Source

| Source | Meaning |
|--------|---------|
| `TLOG_BASED` | The log file was parsed and the actual max-fractional-SEM value was verified |
| `SUMMARY_BASED` | Only p.out parameters were available (iteration count, restart count) |
| `INSUFFICIENT` | Key parameters were missing |

These labels describe **what data was used** for the assessment.  They are
NOT statistical confidence levels.

## Iteration Status

| Status | Meaning |
|--------|---------|
| `MAX_ITERATION_CAP_REACHED` | The engine ran all configured iterations |
| `STOPPED_BEFORE_CAP` | The engine stopped before the iteration limit |
| `UNKNOWN` | Iteration information not available |
| `NOT_APPLICABLE` | No outer Gibbs iteration applies to the selected mode |

Iteration status is **separate** from stability status.  The engine can
hit the iteration cap while still satisfying the stability criterion.

## Example: Phase 5 Real Run

```
max_fractional_sem = 0.00606
threshold          = 0.01
stability_status   = PIGEAN_STABILITY_CRITERION_MET
assessment_source  = TLOG_BASED
iteration_status   = MAX_ITERATION_CAP_REACHED
formal_convergence = NOT_FORMALLY_ASSESSED
```

The engine ran all 500 iterations (hit cap).  At completion, the
max-fractional-SEM was 0.00606, which is below the 0.01 threshold.
The stability criterion was met despite the iteration cap.

Formal MCMC convergence was not assessed because the engine does not
retain individual chain traces.

## Threshold Source

The threshold `0.01` comes from the original PIGEAN engine:

```python
# priors.py line 498
parser.add_option("","--max-frac-sem",type=float,default=0.01)
```

The Phase 5 wrapper uses this same threshold.  It was not changed.

## Audit Reference

The scientific audit that led to this terminology correction is preserved at:

```
outputs/phase5_report/convergence/convergence_audit.md
```

Classification: **LEVEL 2 — ENGINE-SPECIFIC STABILITY CRITERION**
Recommendation: **B — KEEP THE CALCULATION, RENAME THE STATUS**
