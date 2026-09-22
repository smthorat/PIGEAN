# GWAS Input — User Guide

## Overview

The PIGEAN wrapper supports GWAS (Genome-Wide Association Study) summary statistics
as input evidence. Unlike gene-level evidence types (positive controls, Bayes factors,
Z-scores, percentiles, exomes), GWAS data is **SNP-level**: each row represents a
genetic variant with its chromosomal coordinates and association statistics.

The PIGEAN engine maps SNPs to genes using genomic coordinates from reference files,
computing HuGE (Human Gene Expression) scores internally. The wrapper validates your
input, enforces genome-build consistency, and routes the correct reference files.

## Genome Build Requirement

**GWAS mode requires an explicit `--genome-build` flag.** Unlike other analysis types,
which default to hg19, the wrapper will refuse to run GWAS without an explicit build
declaration. This prevents silent coordinate-system mismatches.

```bash
# Required for GWAS
python3 run_pigean.py --analysis gwas --genome-build hg19 --input gwas.tsv --output results

# This will ERROR — no genome build specified
python3 run_pigean.py --analysis gwas --input gwas.tsv --output results
```

### Accepted Build Values

| Input | Normalizes to |
|-------|---------------|
| `hg19`, `GRCh37`, `NCBI37` | `hg19` |
| `hg38`, `GRCh38` | `hg38` |

Case-insensitive. Leading/trailing whitespace stripped.

> **Important:** hg38 reference files are not yet available. Specifying hg38 will
> produce a clear error indicating that the reference files are missing. This will be
> supported in a future release.

## Input File Format

The GWAS input file must be a tab-separated (or whitespace-separated) text file with:

1. **A header row** (required)
2. **At least one data row** (required)
3. Columns for chromosome, position, and at least 2 of {p-value, beta, SE, sample size}

### Standard Column Names (Auto-Detected)

The engine auto-detects standard column names. If your file uses these names, no
column override flags are needed:

| Field | Common names detected |
|-------|----------------------|
| Chromosome | `CHR`, `chr`, `chrom`, `chromosome` |
| Position | `BP`, `pos`, `position`, `bp` |
| P-value | `P`, `p`, `pvalue`, `p_value`, `P-value` |
| Beta/effect | `BETA`, `beta`, `b`, `effect` |
| Standard error | `SE`, `se`, `stderr` |
| Sample size | `N`, `n`, `sample_size` |
| Allele frequency | `FRQ`, `freq`, `maf`, `MAF` |

### Custom Column Names

If your file uses non-standard column names, specify them explicitly:

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input my_gwas.tsv \
    --output results \
    --gwas-chrom-col chromosome \
    --gwas-pos-col base_position \
    --gwas-p-col association_pvalue
```

### All Column Override Flags

| Flag | Description |
|------|-------------|
| `--gwas-chrom-col` | Chromosome column name |
| `--gwas-pos-col` | Position column name |
| `--gwas-p-col` | P-value column name |
| `--gwas-beta-col` | Beta/effect size column name |
| `--gwas-se-col` | Standard error column name |
| `--gwas-n-col` | Per-SNP sample size column name |
| `--gwas-n` | Global sample size (scalar, not a column) |
| `--gwas-freq-col` | Allele frequency column name |
| `--gwas-locus-col` | Compound chr:pos locus column (alternative to separate chrom+pos) |
| `--gwas-filter-col` | Row-filter column name |
| `--gwas-filter-value` | Row-filter match value (only rows matching this value are kept) |

Columns can also be specified as 1-based integer indices (e.g., `--gwas-chrom-col 1`
means "use the first column as chromosome").

## How It Works

1. **Validation**: The wrapper checks that the file exists, has a header, has data rows,
   and that any user-specified column names exist in the header.

2. **Normalization**: The file is copied to the output directory. Blank lines are removed.
   No data transformation is performed — the engine handles all computation.

3. **Reference Resolution**: Based on the genome build, the wrapper resolves coordinate
   reference files (gene locations, TSS positions, exon boundaries).

4. **Engine Execution**: The engine receives the GWAS file via `--gwas-in` along with
   genome-build-matched reference files. It performs:
   - SNP-to-gene mapping using distance-based coordinate lookup
   - Per-SNP: p/beta/SE → Z-score → variant Bayes factor
   - Per-gene: aggregate variant BFs → gene-level log-BF (HuGE score)
   - Gene-set enrichment via Gibbs sampling

## Examples

### Basic GWAS (standard column names)

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_sumstats.tsv \
    --output runs/gwas_basic
```

### GWAS with Custom Columns

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build GRCh37 \
    --input ukbb_gwas.tsv \
    --output runs/ukbb \
    --gwas-chrom-col CHROM \
    --gwas-pos-col POS \
    --gwas-p-col PVAL \
    --gwas-beta-col EFFECT \
    --gwas-se-col STDERR
```

### GWAS with Global Sample Size

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_no_n.tsv \
    --output runs/gwas_n \
    --gwas-n 50000
```

### GWAS with Row Filter

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_multiancestry.tsv \
    --output runs/european \
    --gwas-filter-col ANCESTRY \
    --gwas-filter-value EUR
```

## WDL Usage

```wdl
call rock_pigean.run_pigean {
    input:
        input_file = gwas_sumstats,
        analysis_type = "gwas",
        genome_build = "hg19",
        gwas_chrom_col = "CHR",
        gwas_pos_col = "BP",
}
```

## Troubleshooting

### "GWAS analysis requires --genome-build to be specified explicitly"

You must pass `--genome-build` when using `--analysis gwas`. The wrapper will not
default to hg19 for GWAS to prevent silent coordinate mismatches.

### "Required reference file not found"

The requested genome build's reference files are not installed. Currently only hg19
reference files are available. If you're using hg38 coordinates, you'll need to wait
for hg38 reference file support or lift over your GWAS to hg19 coordinates externally.

### "Specified column 'X' not found in header"

The column name you specified via `--gwas-*-col` doesn't match any column in your
file's header. Check for typos or case mismatches.

## Technical Notes

- The engine's chromosome handling strips "chr" prefixes (e.g., "chr1" → "1").
  Both "chr1" and "1" formats work.
- P-values are clamped to ≥1e-250 by the engine.
- "NA" values in any field are silently skipped.
- If SE is missing but N is provided, SE is derived as 1/sqrt(N).
- Duplicate SNPs at the same position are all kept.
- Gene QC is not applicable for GWAS — the engine maps SNPs to genes internally.
  The wrapper reports this in the QC section.
