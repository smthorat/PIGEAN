# Evidence Input Types — PIGEAN Wrapper

This document describes the gene-level evidence types supported by the PIGEAN wrapper
(Phase 3: Multi-Evidence Input Adapters).

## Overview

The wrapper supports five analysis/evidence types, selectable via `--analysis`:

| Analysis Type | `--analysis` Value | Input Format | Required Columns |
|---|---|---|---|
| Positive Controls | `positive-controls` (default) | One gene per line | — |
| Gene Bayes Factors | `gene-bayes-factor` | TSV with header | `Gene`, `log_bf` |
| Gene Z-Scores | `gene-z-score` | TSV with header | User-specified |
| Gene Percentiles | `gene-percentile` | TSV with header | User-specified |
| Exome Associations | `exome` | TSV with header | Auto-detected |

All types produce the same output files: `gs.out`, `gss.out`, `ggss.out`, `p.out`.

---

## 1. Positive Controls (`--analysis positive-controls`)

**Default mode.** Input is a plain text file with one gene symbol per line.

```bash
python3 run_pigean.py --input gene_list.txt --output results/
```

The wrapper normalizes the gene list (deduplicates, trims whitespace, removes blanks)
and passes it to the engine via `--positive-controls-in`.

**Interpretation:** Positive control genes are INPUT genes, not independently discovered.
Their high priors are expected by design. Additional genes with elevated priors are
prioritized candidates based on gene-set co-membership patterns.

---

## 2. Gene Bayes Factors (`--analysis gene-bayes-factor`)

Input is a TSV file with gene-level natural-log Bayes factors.

```bash
python3 run_pigean.py \
    --analysis gene-bayes-factor \
    --input evidence.tsv \
    --output results/
```

### Input Format

Default columns: `Gene` and `log_bf`. Column names can be overridden:

```bash
python3 run_pigean.py \
    --analysis gene-bayes-factor \
    --input evidence.tsv \
    --gene-column GENE_ID \
    --score-column my_log_bf \
    --output results/
```

### Scale Convention

**Values MUST be natural-log Bayes factors (ln(BF)).**

- Positive values = evidence for association
- Negative values = evidence against association
- Zero = no evidence

The engine reads these values directly without transformation. If your values are
log10(BF), convert them first: `ln(BF) = log10(BF) × ln(10) ≈ log10(BF) × 2.3026`.

### Duplicate Handling

The engine uses **last-value-wins** for duplicate gene entries.

---

## 3. Gene Z-Scores (`--analysis gene-z-score`)

Input is a TSV file with gene-level score values. Both `--gene-column` and
`--score-column` are **required**.

```bash
python3 run_pigean.py \
    --analysis gene-z-score \
    --input scores.tsv \
    --gene-column Gene \
    --score-column Score \
    --output results/
```

### ⚠️ Critical Semantic Note

Despite the name "z-score", **the engine does NOT treat these as statistical
Z-scores (beta/SE)**. Instead, it treats them as **raw log-odds values**:

```
bf = value - mean(values) + background_log_bf
```

The engine:
1. Mean-centers the input values
2. Shifts by a background log-BF term
3. Uses the result as log-BF evidence

This is appropriate for gene-level scores on a log-odds-like scale.

### Duplicate Handling

The engine uses **last-value-wins** for duplicate gene entries.

---

## 4. Gene Percentiles (`--analysis gene-percentile`)

Input is a TSV file with gene-level rank/percentile scores. Both `--gene-column`
and `--score-column` are **required**.

```bash
python3 run_pigean.py \
    --analysis gene-percentile \
    --input ranks.tsv \
    --gene-column Gene \
    --score-column rank_score \
    --output results/
```

### Direction Flag

By default, **lower values = stronger evidence** (e.g., p-values, ranks).
If higher values indicate stronger evidence, add `--higher-is-better`:

```bash
python3 run_pigean.py \
    --analysis gene-percentile \
    --input ranks.tsv \
    --gene-column Gene \
    --score-column enrichment_score \
    --higher-is-better \
    --output results/
```

### Internal Conversion

The engine converts percentiles to log-BFs via:
1. Rank → quantile
2. Quantile → inverse normal CDF
3. Result used as log-BF evidence

### Duplicate Handling

The engine uses **last-value-wins** for duplicate gene entries.

---

## 5. Exome Associations (`--analysis exome`)

Input is a TSV file with gene-level exome/rare-variant association statistics.

```bash
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --output results/
```

### Column Auto-Detection

The engine auto-detects columns by searching for standard headers:
- **Gene:** `Gene`, `gene`, `GENE`, `gene_id`, …
- **P-value:** `p`, `P`, `pval`, `Pvalue`, …
- **Beta:** `beta`, `Beta`, `BETA`, `b`, …
- **SE:** `se`, `SE`, `Se`, `std_error`, …
- **N:** `n`, `N`, `sample_size`, …

You can override column names:

```bash
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --exomes-gene-col GENE_SYMBOL \
    --exomes-p-col PVALUE \
    --exomes-n 50000 \
    --output results/
```

### Required Statistics

The engine requires **at least 2 of 3** statistics per gene: p-value, beta, SE.
From these it computes: Z → HuGE score (log-BF).

### Available Overrides

| Flag | Description |
|---|---|
| `--exomes-gene-col` | Gene column name |
| `--exomes-p-col` | P-value column name |
| `--exomes-beta-col` | Beta/effect size column name |
| `--exomes-se-col` | Standard error column name |
| `--exomes-n-col` | Sample size column name |
| `--exomes-n` | Global sample size (scalar) |

### Duplicate Handling

The engine uses **first-occurrence-wins** for duplicate genes (unlike BF/Z-score modes).

---

## Validation

All evidence types undergo two validation stages:

1. **File validation** — checks file exists, is not empty, has required columns
2. **Gene QC** — checks gene symbols against the gene map

Results are recorded in:
- `input_gene_qc.tsv` (positive-controls mode)
- `input_evidence_qc.tsv` (all other modes)
- `run_manifest.json` (full provenance)

### Validation Status

| Status | Meaning |
|---|---|
| `PASS` | All checks passed |
| `PASS_WITH_WARNINGS` | Non-fatal issues (e.g., unrecognized genes) |
| `FAIL` | Critical issues — engine will NOT run |

---

## WDL Usage

All evidence parameters are exposed in the WDL workflow:

```json
{
    "rock_pigean.input_file": "gs://bucket/evidence.tsv",
    "rock_pigean.analysis_type": "gene-bayes-factor",
    "rock_pigean.gene_column": "Gene",
    "rock_pigean.score_column": "log_bf"
}
```

See `wdl/rock_pigean.wdl` for the full input schema.

---

## Not Supported

The following evidence types are NOT exposed through the wrapper:

- **GWAS summary statistics** — requires SNP-level processing (Phase 4)
- **Probability input** — requires sigmoid transformation (not traced)
- **Z-score logistic mode** — requires `--gene-zs-logistic` (advanced)
- **Credible sets** — requires additional file format handling

These may be added in future phases.
