# PIGEAN Pipeline — Results Guide

This document describes the four output files produced by the PIGEAN (Probabilistic Inference of Gene ENrichment from Association Networks) pipeline and how to interpret them.

---

## Quick Overview

| Output File | Description | Rows | Key Use |
|---|---|---|---|
| `gs.*.out` | Gene-level statistics | One per gene (~42,216) | Which genes are most likely relevant? |
| `gss.*.out` | Gene-set-level statistics | One per gene set (~41,627) | Which pathways are enriched? |
| `ggss.*.out` | Gene × gene-set pair statistics | One per gene–pathway pair (~283K) | Why does a gene have a high prior? |
| `p.*.out` | Model hyperparameters | ~36 rows | What settings/parameters were learned? |

---

## 1. Gene Stats — `gs.*.out`

**What it is:** Per-gene scores estimating how likely each gene is to be relevant to your input gene list, based on shared pathway membership.

### Column Definitions

| Column | Description |
|---|---|
| `Gene` | Gene symbol |
| `prior` | **Pathway-based prior probability** — how much the gene sets alone predict this gene's relevance. Higher = more pathway support. Range: negative (no support) to ~2+ (strong support) |
| `prior_adj` | Prior after adjustment for genomic covariates (gene size, location bias) |
| `combined` | **Combined score** = `prior` + `log_bf`. For input genes, this is dominated by `log_bf`. For non-input genes, this ≈ `prior` |
| `combined_adj` | Combined score after adjustment |
| `combined_D` | **Posterior probability** that the gene is relevant (0–1 scale). This is the most interpretable single number |
| `positive_control` | Log Bayes factor from being in the input gene list. `5.89` = input gene; `NA` = not in input list |
| `log_bf` | Log Bayes factor from the evidence (same as `positive_control` for this pipeline mode) |
| `N` | Number of gene sets this gene belongs to |
| `Chrom` | Chromosome |
| `Start` | Genomic start position (hg19) |
| `End` | Genomic end position (hg19) |

### Snapshot

```
Gene       prior   prior_adj  combined  combined_adj  combined_D  positive_control  log_bf  N     Chrom  Start       End
LEP        1.89    1.92       7.78      7.79          0.986       5.89              5.89    930   7      127881241   127897682
ADIPOQ     1.36    1.38       7.25      7.26          0.979       5.89              5.89    650   3      186560463   186576252
BRD2       0.803   0.823      6.69      6.7           0.969       5.89              5.89    432   6      32936437    32949282
SOX2       0.738   0.759      6.63      6.64          0.967       5.89              5.89    499   3      181429712   181432224
LITAF      0.223   0.243      6.11      6.12          0.957       5.89              5.89    453   16     11641578    11681322
LECT2      0.197   0.212      6.09      6.1           0.955       5.89              5.89    114   5      135282600   135290723
SLCO1B1    0.175   0.19       6.06      6.08          0.955       5.89              5.89    119   12     21284128    21392730
RMI2       0.0481  0.0641     5.94      5.95          0.951       5.89              5.89    164   16     11439311    11445617
CD300LG    0.0383  0.0528     5.93      5.94          0.951       5.89              5.89    66    17     41924145    41940997
PLINK1     -0.0104 0.00309    5.88      5.89          0.95        5.89              5.89    0     NA     NA          NA
LEPR       0.988   1.01       0.991     1             0.174       NA                0       677   1      65886335    66103176
PPARG      0.64    0.669      0.64      0.65          0.124       NA                0       1003  3      12329349    12475855
APOE       0.432   0.465      0.434     0.443         0.097       NA                0       1238  19     45409039    45412650
```

### How to Interpret

- **Input genes** (positive_control = 5.89): All have high combined scores (5.88–7.78) because the log Bayes factor dominates. The `prior` column tells you how well the pathways independently support each gene. LEP (prior 1.89) has the most pathway evidence; PLINK1 (prior -0.01) has none.
- **Non-input genes** (positive_control = NA): The `prior` column is the key metric. These are "discovered" genes — genes the model predicts should be related based on pathway overlap. LEPR (0.988), PPARG (0.64), and APOE (0.43) are the top discoveries.
- **`combined_D`** is the easiest number to use: it's the posterior probability (0–1) that the gene is relevant. Values > 0.5 suggest strong evidence; values near the background (~0.05) suggest no evidence.
- **`N`** tells you how many gene sets the gene belongs to. Genes with very low N (like PLINK1 with 0) can't get pathway support because they aren't in the reference databases.

---

## 2. Gene Set Stats — `gss.*.out`

**What it is:** Per-pathway scores estimating how associated each gene set (pathway) is with your input gene list.

### Key Column Definitions

| Column | Description |
|---|---|
| `Gene_Set` | Gene set / pathway name |
| `label` | Source database (`gene_set_list_mouse_2024` or `gene_set_list_msigdb_nohp`) |
| `N` | Number of genes in this gene set |
| `scale` | Scaling factor based on gene set size |
| `beta_tilde` | **Marginal effect size** (univariate association with your gene list) |
| `P` | P-value for the marginal association |
| `Z` | Z-score for the marginal association |
| `SE` | Standard error of the beta |
| `beta` | **Conditional effect size** — the gene set's effect after accounting for all other gene sets. This is the most important score |
| `beta_uncorrected` | Beta before correction for correlated gene sets |
| `avg_postp` | Average posterior probability across member genes |
| `*_orig` columns | Same statistics but from the original (pre-Gibbs) estimation |
| `p_used` | Learned prior probability hyperparameter used for this gene set's batch |
| `sigma2_used` | Learned variance hyperparameter used for this gene set's batch |

### Snapshot

```
Gene_Set                                          label                    N    beta_tilde  P          Z     beta      beta_uncorrected
mp_increased_pancreatic_islet_number              gene_set_list_mouse_2024 12   1.51        0.0289     2.18  0.0224    0.0589
mp_increased_interscapular_fat_pad_weight         gene_set_list_mouse_2024 15   1.4         0.0277     2.2   0.0191    0.0989
mp_increased_fat_cell_size                        gene_set_list_mouse_2024 50   0.661       0.124      1.54  0.0145    0.0615
mp_abnormal_gluconeogenesis                       gene_set_list_mouse_2024 55   0.639       0.123      1.54  0.0144    0.0607
WP_LEPTIN_AND_ADIPONECTIN                         gene_set_list_msigdb_nohp 10  1.75        0.0133     2.48  0.00938   0.0742
mp_increased_food_intake                          gene_set_list_mouse_2024 145  0.44        0.118      1.56  0.01      0.0641
```

### How to Interpret

- **`beta`** (conditional): The most reliable column. It represents the independent contribution of each pathway after accounting for overlap with other pathways. Higher values mean the pathway independently explains membership in your gene list.
- **`beta_tilde`** (marginal): The raw association before conditioning. Useful for initial filtering but can be inflated by correlated gene sets.
- **`P` and `Z`**: Significance of the marginal association. Low P / high Z = pathway is enriched in your gene list. Note: these are from the marginal model, not the conditional.
- **`beta_uncorrected` vs `beta`**: The difference shows how much overlap correction matters. Large gaps mean this pathway's signal is partially explained by other pathways.
- **`label`**: The source database. `gene_set_list_mouse_2024` = mouse phenotype gene sets; `gene_set_list_msigdb_nohp` = MSigDB (human canonical pathways, GO terms, curated sets).

---

## 3. Gene–Gene-Set Stats — `ggss.*.out`

**What it is:** The "explanation file" — for each gene, which gene sets contributed to its prior, and by how much.

### Column Definitions

| Column | Description |
|---|---|
| `Gene` | Gene symbol |
| `prior` | Gene's overall prior (same as in gs.out) |
| `combined` | Gene's combined score (same as in gs.out) |
| `log_bf` | Gene's log Bayes factor |
| `gene_set` | The pathway contributing to this gene's prior |
| `beta` | The gene set's conditional effect size (how much this pathway contributes) |
| `weight` | Gene's membership weight in this gene set (1 = full member) |

### Snapshot — LEP (high prior gene)

```
Gene  prior  combined  log_bf  gene_set                                       beta      weight
LEP   1.89   7.78      5.89    mp_increased_pancreatic_islet_number           0.0224    1
LEP   1.89   7.78      5.89    mp_increased_interscapular_fat_pad_weight      0.0191    1
LEP   1.89   7.78      5.89    mp_increased_fat_cell_size                     0.0145    1
LEP   1.89   7.78      5.89    mp_abnormal_gluconeogenesis                    0.0144    1
LEP   1.89   7.78      5.89    GOBP_PROSTAGLANDIN_TRANSPORT                   0.0143    1
LEP   1.89   7.78      5.89    mp_decreased_carbon_dioxide_production         0.0135    1
LEP   1.89   7.78      5.89    mp_decreased_fatty_acid_beta-oxidation         0.0105    1
LEP   1.89   7.78      5.89    WP_LEPTIN_AND_ADIPONECTIN                      0.00938   1
```

### Snapshot — LECT2 (moderate prior gene)

```
Gene   prior  combined  log_bf  gene_set                                      beta      weight
LECT2  0.197  6.09      5.89    mp_abnormal_leukocyte_physiology              0.00646   1
LECT2  0.197  6.09      5.89    mp_increased_hepatocyte_apoptosis             0.00568   1
LECT2  0.197  6.09      5.89    mp_liver_inflammation                         0.00558   1
LECT2  0.197  6.09      5.89    mp_increased_NK_T_cell_number                 0.00423   1
LECT2  0.197  6.09      5.89    mp_abnormal_circulating_enzyme_level          0.00305   1
LECT2  0.197  6.09      5.89    mp_hepatocellular_carcinoma                   0.00263   1
```

### How to Interpret

- Each row shows one gene–pathway link. The `beta` value is that pathway's contribution to the gene's overall `prior`.
- **A gene's prior ≈ sum of (beta × weight) across all its gene-set entries.** This lets you trace exactly why a gene scored high or low.
- LEP's high prior (1.89) comes from many metabolic pathways each contributing a small amount — fat cell size, gluconeogenesis, adipose tissue, food intake, leptin signaling.
- LECT2's moderate prior (0.197) comes from fewer, smaller contributions — hepatic and immune pathways.
- Use this file to answer: *"Why did gene X score high/low?"*

---

## 4. Parameters — `p.*.out`

**What it is:** The model hyperparameters — both fixed settings and values learned during Gibbs sampling.

### Snapshot

```
Parameter                        Version  Value
p                                1        0.00356
p                                2        0.00187
sigma2                           1        6.41e-08
sigma2                           2        3.36e-08
sigma2_cond                      1        1.80e-05
sigma2_cond                      2        1.80e-05
num_gene_sets_read               1        2094
num_genes_read                   1        42216
num_chains                       1        10
num_gibbs_iter                   1        499
sparse_solution                  1        True
sparse_frac                      1        0.01
max_num_gene_sets                1        5000
filter_gene_set_p                1        0.01
adjust_priors                    1        True
```

### Key Parameters

| Parameter | Meaning |
|---|---|
| `p` | **Learned proportion of gene sets with non-zero effects.** Version 1 and 2 correspond to the two gene set input files. ~0.2–0.4% of gene sets have effects |
| `sigma2` | **Learned effect size variance.** Controls how large gene set betas can be |
| `sigma2_cond` | Conditional variance (fixed based on top gene set prior) |
| `num_gene_sets_read` | Total gene sets retained after filtering (from ~41K input) |
| `num_genes_read` | Total genes in the analysis |
| `num_chains` | Number of parallel MCMC chains (10) |
| `num_gibbs_iter` | Number of Gibbs sampling iterations completed (499) |
| `sparse_solution` | Whether sparse beta estimation was used |
| `filter_gene_set_p` | P-value threshold for gene set inclusion |

---

## Practical Tips

### Finding the top genes predicted by the model

```bash
# Top 20 non-input genes by prior (discovered genes)
awk -F'\t' 'NR>1 && $7=="NA"' gs.*.out | sort -t$'\t' -k2 -rn | head -20

# All input genes ranked by prior
awk -F'\t' 'NR>1 && $7!="NA"' gs.*.out | sort -t$'\t' -k2 -rn
```

### Finding the top enriched pathways

```bash
# Top 20 pathways by conditional beta
awk -F'\t' 'NR>1' gss.*.out | sort -t$'\t' -k10 -rn | head -20

# Top pathways by marginal significance
awk -F'\t' 'NR>1' gss.*.out | sort -t$'\t' -k7 -g | head -20
```

### Explaining a specific gene's score

```bash
# What pathways drive LEP's high prior?
grep "^LEP	" ggss.*.out | sort -t$'\t' -k6 -rn | head -10

# What pathways drive LECT2's prior?
grep "^LECT2	" ggss.*.out | sort -t$'\t' -k6 -rn | head -10
```

### Filtering for high-confidence results

```bash
# Genes with posterior probability > 0.5 (strong evidence)
awk -F'\t' 'NR>1 && $6 > 0.5' gs.*.out

# Non-input genes with prior > 0.1 (pathway-predicted novel genes)
awk -F'\t' 'NR>1 && $7=="NA" && $2 > 0.1' gs.*.out | sort -t$'\t' -k2 -rn
```

---

## Input Gene List Used

The gene list provided as input (`--positive-controls-in`):

```
RMI2
LITAF
SLCO1B1
BRD2
SOX2
LEP
PLINK1
ADIPOQ
LECT2
CD300LG
```

## Reference Data (hg19/GRCh37)

All genomic coordinates use **hg19/GRCh37**. The pipeline uses two gene set databases:

- **Mouse phenotype gene sets** (`gene_set_list_mouse_2024.txt`) — gene sets derived from mouse phenotype annotations
- **MSigDB gene sets** (`gene_set_list_msigdb_nohp.txt`) — Molecular Signatures Database (canonical pathways, GO terms, curated gene sets), excluding HP (Human Phenotype) terms
