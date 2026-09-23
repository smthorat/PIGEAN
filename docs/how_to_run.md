# How to Run PIGEAN

This guide is for researchers running the repository from a fresh clone. It
shows how to select an evidence type, provide either built-in or custom gene
sets, run locally or with Docker, and verify the result.

## What the One-Command Workflow Does

One `run_pigean.py` command runs the complete wrapper workflow for one selected
analysis:

1. validates the input and requested options;
2. preserves the original input and creates an engine-ready normalized copy;
3. performs input and gene-mapping quality control;
4. resolves gene-set and genome-reference files;
5. runs the PIGEAN scientific engine;
6. parses the engine outputs and assesses run stability;
7. writes researcher-facing reports and a provenance manifest.

The command does not guess the scientific meaning of a file. The researcher
must select the correct `--analysis` and, when column names are nonstandard,
identify those columns explicitly.

## Current Capability

| Evidence supplied by the researcher | `--analysis` | `standard` | `naive-priors` |
|---|---|---:|---:|
| One gene symbol per line | `positive-controls` | Supported | Supported |
| Gene-level natural-log Bayes factors | `gene-bayes-factor` | Supported | Not verified; rejected |
| Gene-level exome/rare-variant statistics | `exome` | Supported | Not verified; rejected |
| SNP-level GWAS summary statistics | `gwas` | Supported | Not verified; rejected |
| Gene-level Z scores | `gene-z-score` | Engine-blocked; rejected | Engine-blocked; rejected |
| Gene-level percentiles | `gene-percentile` | Engine-blocked; rejected | Engine-blocked; rejected |

`standard` is the default and should be used unless the study specifically
requires the verified positive-control-only `naive-priors` estimator. See
[Advanced Modes](advanced_modes.md) for the scientific boundary between them.

The command-line help still lists Z-score and percentile adapters because their
input contracts are retained, but the wrapper deliberately rejects them before
engine execution. Do not use those two analysis values for production runs.

## 1. Get the Repository

```bash
git clone https://github.com/smthorat/PIGEAN.git
cd PIGEAN
```

Run all local commands from the repository root. Relative paths in the examples
are resolved from that directory.

## 2. Choose a Runtime

### Option A: Docker (recommended for reproducibility)

Docker packages Python, NumPy, SciPy, the wrapper, the frozen engine, and the
included hg19 reference data into one image.

```bash
docker build -t rock-pigean:6.0.0 .
```

Create an output-parent directory on the host, then mount the input directory
read-only and the output directory read-write:

```bash
mkdir -p /absolute/path/to/output-parent

docker run --rm \
  -v /absolute/path/to/data:/inputs:ro \
  -v /absolute/path/to/output-parent:/outputs \
  rock-pigean:6.0.0 \
  --analysis positive-controls \
  --input /inputs/my_genes.txt \
  --output /outputs/my_positive_controls
```

Replace `/absolute/path/to/data` and `/absolute/path/to/output-parent` with
real absolute paths. Docker Desktop users must allow Docker to access the
mounted host directories.

The arguments after the image name are the same arguments used in every local
example below. Files such as custom gene sets and backgrounds must also be
inside a mounted directory and referenced by their container paths.

### Option B: Local Python

Local execution requires Python 3, NumPy, and SciPy:

```bash
python3 -c "import numpy, scipy; print('Dependencies available')"
python3 run_pigean.py --help
```

If the import fails, install NumPy and SciPy in the Python environment selected
for the analysis. The repository already contains the engine and the included
reference data; they do not need to be downloaded separately.

## 3. Run the Scenario That Matches the Input

Use a new output directory for every run. The wrapper refuses to write into a
nonempty output directory unless `--overwrite` is supplied.

### Scenario A: Positive-control gene list

The file is plain text with one gene symbol per line and no header:

```text
LEP
ADIPOQ
BRD2
```

Run the default, validated standard workflow:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --mode standard \
  --input data/my_genes.txt \
  --output results/runs/my_positive_controls
```

Blank lines and duplicate genes are removed from the normalized copy. The
original file is preserved in the output directory.

To use the alternative verified estimator for positive controls:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --mode naive-priors \
  --input data/my_genes.txt \
  --output results/runs/my_naive_priors
```

`naive-priors` does not run the outer gene-prior Gibbs loop and does not emit
`combined_D`. Its report ranks genes using `combined`; its outer-loop stability
and formal MCMC-convergence fields are marked `NOT_APPLICABLE`.

### Scenario B: Gene-level Bayes factors

Provide a whitespace- or tab-delimited file with a header. The default columns
are `Gene` and `log_bf`:

```text
Gene	log_bf
LEP	2.5
ADIPOQ	1.8
BRD2	-0.4
```

The score must be the natural logarithm of the Bayes factor, not log10(BF).

```bash
python3 run_pigean.py \
  --analysis gene-bayes-factor \
  --mode standard \
  --input data/gene_bayes_factors.tsv \
  --output results/runs/my_bayes_factors
```

For different column names, map them explicitly:

```bash
python3 run_pigean.py \
  --analysis gene-bayes-factor \
  --input data/gene_bayes_factors.tsv \
  --gene-column gene_symbol \
  --score-column ln_bf \
  --output results/runs/my_bayes_factors
```

Nonfinite or nonnumeric scores are removed with QC warnings. If a gene occurs
more than once, the engine uses its last occurrence.

### Scenario C: Exome or rare-variant association results

Provide a whitespace- or tab-delimited file with a header, a gene column, and
at least two of p-value, beta, and standard-error columns. Common column names
are detected automatically. A typical input is:

```text
Gene	P	BETA	SE	N
LEP	0.0001	0.40	0.10	50000
ADIPOQ	0.0100	-0.25	0.09	50000
```

```bash
python3 run_pigean.py \
  --analysis exome \
  --mode standard \
  --input data/exome_results.tsv \
  --output results/runs/my_exome_run
```

For nonstandard headers, provide the relevant mappings:

```bash
python3 run_pigean.py \
  --analysis exome \
  --input data/exome_results.tsv \
  --exomes-gene-col SYMBOL \
  --exomes-p-col PVALUE \
  --exomes-beta-col EFFECT \
  --exomes-se-col STDERR \
  --exomes-n-col SAMPLE_SIZE \
  --output results/runs/my_exome_run
```

If sample size is not a column, a single study-wide value can be supplied with
`--exomes-n 50000`. The engine derives the missing member of p-value, beta, and
SE when enough information is present. For duplicate genes, the first
occurrence is retained by the engine.

### Scenario D: GWAS summary statistics

GWAS input is variant-level, must have a header, and must use hg19/GRCh37
coordinates with the currently included reference files. Provide either
chromosome and position columns or a compound locus column, plus at least two
of p-value, beta, standard error, and sample size (a global sample size also
counts).

A typical input is:

```text
CHR	BP	P	BETA	SE	N
1	123456	0.00001	0.15	0.03	100000
1	234567	0.02000	-0.08	0.04	100000
```

The genome build is mandatory in GWAS mode:

```bash
python3 run_pigean.py \
  --analysis gwas \
  --mode standard \
  --genome-build hg19 \
  --input data/gwas_sumstats.tsv \
  --output results/runs/my_gwas_run
```

For nonstandard headers:

```bash
python3 run_pigean.py \
  --analysis gwas \
  --genome-build GRCh37 \
  --input data/gwas_sumstats.tsv \
  --gwas-chrom-col CHROM \
  --gwas-pos-col POS \
  --gwas-p-col PVAL \
  --gwas-beta-col EFFECT \
  --gwas-se-col STDERR \
  --gwas-n 100000 \
  --output results/runs/my_gwas_run
```

GWAS column selectors may be header names or 1-based column indices. Both `1`
and `chr1` chromosome values are accepted. Although the CLI recognizes hg38 and
GRCh38 aliases, hg38 reference files are not included, so an hg38 run will stop
with a missing-reference error. See [GWAS Input](gwas_input.md) for filtering,
locus-column, allele-frequency, and other GWAS-specific options.

## 4. Select Gene Sets

The built-in profiles are:

| Value | Gene sets used |
|---|---|
| `default` | Mouse annotations plus MSigDB annotations |
| `mouse-only` | Mouse annotations only |
| `msigdb-only` | MSigDB annotations only |
| `custom` | Researcher-supplied files; behavior depends on replace/supplement |

Select a built-in profile by adding one flag to any supported run:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --input data/my_genes.txt \
  --gene-sets msigdb-only \
  --output results/runs/msigdb_only
```

### Add custom gene sets to a built-in profile

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --input data/my_genes.txt \
  --gene-sets default \
  --custom-gene-set-files data/lab_pathways.gmt \
  --custom-gene-set-format gmt \
  --custom-gene-set-action supplement \
  --output results/runs/default_plus_lab
```

### Use only custom gene sets

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --input data/my_genes.txt \
  --gene-sets custom \
  --custom-gene-set-files data/lab_pathways.gmt \
  --custom-gene-set-format gmt \
  --custom-gene-set-action replace \
  --output results/runs/custom_only
```

Multiple custom files may follow `--custom-gene-set-files`. Every custom-file
run must specify both its format and its action.

Supported custom formats:

```text
# GMT: set name, description, then genes
MY_PATHWAY	description	LEP	ADIPOQ	BRD2

# Engine format: set name followed immediately by genes
MY_PATHWAY	LEP	ADIPOQ	BRD2
```

The wrapper preserves each original file and converts GMT to engine format in
the run's `normalized/` directory.

## 5. Optional Positive-Control Background

A background is a second one-gene-per-line file representing the set from
which positive controls could have been selected. Use it with a
positive-control analysis:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --input data/my_genes.txt \
  --background data/assayed_genes.txt \
  --output results/runs/positive_controls_with_background
```

The background should contain all input positive-control genes. A background
with no recognized genes fails validation; missing input genes are reported as
a QC warning. The original and normalized background files are preserved.

## 6. Reuse Settings with a JSON Config

Settings can be stored in JSON while file paths remain explicit on the command
line. For example, `config/bayes_factor.json` could contain:

```json
{
  "analysis": "gene-bayes-factor",
  "mode": "standard",
  "gene_sets": "default",
  "gene_column": "gene_symbol",
  "score_column": "ln_bf"
}
```

Run it with:

```bash
python3 run_pigean.py \
  --config config/bayes_factor.json \
  --input data/gene_bayes_factors.tsv \
  --output results/runs/configured_bayes_factors
```

Explicit command-line settings override values from the JSON file. Keep the
primary input path, output path, background path, and custom gene-set file paths
on the command line rather than in JSON.

## 7. Optional Standard-Mode Stability Trace

The standard workflow always writes `convergence.json`. To additionally retain
the engine's convergence trace, add:

```bash
--enable-convergence-trace
```

This option is incompatible with `--mode naive-priors`, because that mode does
not run the outer Gibbs pathway being traced. Stability terminology and limits
are documented in [Convergence and Stability](convergence.md).

## 8. Confirm That a Run Succeeded

A successful run exits with code 0 and prints the result location. Inspect these
files first:

| File | Purpose |
|---|---|
| `report.html` | Primary self-contained researcher report |
| `report.txt` | Plain-text version of the report |
| `run_manifest.json` | Validation, execution, mode, provenance, and checksums |
| `convergence.json` | PIGEAN stability assessment and applicability statements |
| `input_gene_qc.tsv` or `input_evidence_qc.tsv` | Input mapping and QC details |
| `resolved_config.json` | Exact effective settings after defaults/config/CLI merge |
| `priors_command.txt` | Exact engine invocation assembled by the wrapper |
| `pigean_run.log` | Engine stdout and stderr |
| `gs.out`, `gss.out`, `ggss.out`, `p.out` | Raw engine result tables |

The standard workflow may report whether the engine-specific PIGEAN stability
criterion was met, but formal multi-chain MCMC convergence is
`NOT_FORMALLY_ASSESSED`. These are intentionally different claims. See the
[Output Contract](output_contract.md) before interpreting raw result columns.

Because the engine uses unseeded NumPy sampling, exact numerical results can
vary between repeated runs. Preserve `run_manifest.json`, `resolved_config.json`,
and the input copies when recording an analysis.

## 9. Terra or Cromwell

The parameterized workflow is [rock_pigean.wdl](../workflows/rock_pigean.wdl).
It exposes the same analysis, mode, gene-set, background, exome, GWAS, and
stability-trace settings. Its default runtime request is 4 CPUs, 8 GB memory,
and a 50 GB local disk. Ensure the workflow's Docker image is available to the
execution environment before launching it.

## Troubleshooting

### The output directory already exists

Choose a new directory for each scientific run. If replacement is intentional,
rerun with `--overwrite` after confirming the path is correct.

### A required column is not found

Check spelling and capitalization in the header, then pass the appropriate
column flag. Bayes-factor inputs use `--gene-column` and `--score-column`;
exome and GWAS inputs use their `--exomes-*` and `--gwas-*` flags.

### Input genes are not recognized

Review the QC file and confirm that the input identifiers match the gene
symbols used by the included gene map. The wrapper records recognized and
unrecognized values rather than silently treating them as equivalent.

### GWAS reports missing reference files

Confirm that the input coordinates are hg19/GRCh37 and pass
`--genome-build hg19` or `--genome-build GRCh37`. The repository does not
currently include hg38 references.

### Docker cannot read an input file

Use absolute host paths in `-v` mounts, ensure Docker Desktop can access those
directories, and use the mounted container path (for example `/inputs/file.tsv`)
in `--input`.

### The analysis/mode combination is rejected

Use `standard` for positive controls, gene Bayes factors, exome, or GWAS. Use
`naive-priors` only for positive controls. Z-score and percentile evidence are
currently engine-blocked in both modes.

### More command details are needed

```bash
python3 run_pigean.py --help
```

For Docker:

```bash
docker run --rm rock-pigean:6.0.0 --help
```
