# Naive-priors Mode

`naive-priors` is a verified positive-control workflow backed by the engine's
canonical `naive_priors` subcommand. It is not a deterministic or purely
algebraic shortcut: the engine first estimates gene-set effects with its inner
sampling routine, then bypasses only the outer gene-prior Gibbs update.

## Calculation

For each gene, the engine computes:

```text
prior = X_orig · (beta / scale_factor)
prior = prior - mean(prior across present and missing genes)
combined = prior + log_bf
```

With the wrapper's locked settings, adjustment is not applied and
`combined_D` is not produced.

## Inputs and compatibility

- Evidence: positive controls only
- Gene sets: default, mouse-only, msigdb-only, or validated custom profiles
- Background: same normalized positive-control background mechanism as the
  standard wrapper path
- References: gene map plus the selected genome-build locations/TSS/exons and
  selected gene-set files
- Advanced parameters: none

Bayes-factor, exome, and GWAS combinations remain NOT_VERIFIED. Z-score and
percentile combinations retain their existing ENGINE_BLOCKED status.

## Outputs and interpretation

The mode writes `gs.out`, `gss.out`, `ggss.out`, and `p.out`. Reports rank genes
by `combined`, clearly label the prior estimator, and do not display or infer
`combined_D`. The standard and naive priors are not interchangeable.

The outer-Gibbs PIGEAN stability criterion and formal MCMC convergence are both
`NOT_APPLICABLE`. The inner gene-set sampler remains stochastic, which is
reported as a separate limitation.

## Example

```bash
python3 run_pigean.py \
  --mode naive-priors \
  --analysis positive-controls \
  --input examples/gene_list \
  --output results/runs/example_naive
```
