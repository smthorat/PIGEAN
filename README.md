# rock-pigean

PIGEAN (Probabilistic Inference of Gene ENrichment) pipeline with a user-facing wrapper for safer, auditable execution.

## Quick Start

```bash
python3 run_pigean.py \
    --analysis positive-controls \
    --input examples/gene_list \
    --output runs/example_run
```

## Architecture

- **`engine/priors.py`** — the frozen scientific engine (unmodified PIGEAN Gibbs sampler)
- **`run_pigean.py`** — user-facing CLI wrapper
- **`pigean/`** — wrapper modules (config, validation, QC, provenance, convergence, interpretation)

The wrapper adds usability, input QC, convergence assessment, interpretation, and provenance around the existing PIGEAN engine without changing its scientific behavior.

See `docs/architecture.md` for full system design and `docs/repository_structure.md` for directory layout.

## Output Structure

A wrapper run produces:

```
<output_dir>/
├── input/
│   └── positive_controls.txt      # preserved original input
├── normalized/
│   └── positive_controls.txt      # normalized (stripped, deduped)
├── resolved_config.json           # effective scientific settings
├── input_gene_qc.tsv             # per-gene QC results
├── priors_command.txt            # exact engine invocation
├── pigean_run.log                # engine stdout/stderr
├── gs.out                        # gene stats
├── gss.out                       # gene set stats
├── ggss.out                      # gene-gene set stats
├── p.out                         # model parameters
├── convergence.json              # stability assessment
├── run_manifest.json             # full provenance
├── report.txt                    # extended text report
└── report.html                   # self-contained HTML report
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | (required) | Gene list file, one gene per line |
| `--output` | (required) | Output directory |
| `--analysis` | `positive-controls` | Analysis type |
| `--gene-sets` | `default` | Gene set profile |
| `--genome-build` | `hg19` | Genome build |
| `--preset` | `standard` | Parameter preset |
| `--config` | none | Optional JSON config for parameter overrides |
| `--overwrite` | false | Overwrite existing output directory |

## Phase 5: Stability Assessment

Phase 5 adds interpretation reports and PIGEAN stability diagnostics:

- **PIGEAN stability criterion**: Uses the engine's built-in max-fractional-SEM
  check (`--max-frac-sem = 0.01`).  Status is `PIGEAN_STABILITY_CRITERION_MET`
  when the criterion is satisfied.
- **Formal MCMC convergence**: Always `NOT_FORMALLY_ASSESSED` — the frozen
  engine does not retain independent chain traces for R-hat/ESS.
- **Terminology correction**: Following a scientific audit, all statuses use
  precise terminology that distinguishes engine-specific stability from formal
  MCMC convergence.  See `docs/convergence.md` for details.

Additional outputs:
- `convergence.json` — machine-readable stability assessment
- `report.html` — self-contained HTML report with stability section
- `report.txt` — extended text report with top genes, gene sets, metric guide

## Golden Regression Test

```bash
bash tests/golden/run_golden_test.sh              # fresh engine run
bash tests/golden/run_golden_test.sh <output_dir>  # validate existing outputs
```

## WDL (Terra/Cromwell)

See `workflows/rock_pigean.wdl` for the parameterized workflow that calls the wrapper.

## Unit Tests

```bash
python3 -m pytest tests/ -v
```
