# Architecture

## Overview

rock-pigean is a PIGEAN (Probabilistic Inference of Gene ENrichment) pipeline
with a user-facing wrapper for safer, auditable execution. The architecture
follows a strict separation:

- **Frozen engine** (`engine/priors.py`) — the unmodified PIGEAN Gibbs sampler
  (16,567 lines). This file is never changed by the wrapper. SHA-256 is
  recorded in every run manifest.
- **Wrapper** (`run_pigean.py` + `pigean/` package) — adds input validation,
  configuration, QC, provenance tracking, convergence assessment, and
  interpretation reporting.

## Path Resolution

All runtime paths flow through a single function:

```python
# pigean/references.py
def resolve_base_dir():
    if os.path.isfile("/app/engine/priors.py"):  # Docker container
        return "/app"
    # Local development: repo root (parent of pigean/ package)
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
```

From `base_dir`, all paths are resolved:
- Engine: `os.path.join(base_dir, "engine", "priors.py")`
- Reference data: `os.path.join(base_dir, "data", filename)`
- Gene sets: `os.path.join(base_dir, "data", gene_set_file)`

## Data Flow

```
CLI arguments
    │
    ▼
Configuration Resolution (pigean/config.py)
    │ defaults + presets + user overrides
    ▼
Input Validation (pigean/validation.py)
    │ file check → normalize → gene QC
    ▼
Evidence Adapter (pigean/adapters/*.py)
    │ format-specific validation + engine arg building
    ▼
Reference Resolution (pigean/references.py)
    │ genome build → reference paths + gene sets
    ▼
Engine Command Building (pigean/engine.py)
    │ assemble full priors.py command
    ▼
Engine Execution (pigean/engine.py)
    │ subprocess → pigean_run.log
    ▼
Output Parsing (pigean/parsers.py)
    │ gs.out, gss.out, ggss.out, p.out → dicts
    ▼
Convergence Assessment (pigean/convergence.py)
    │ p.out params + log parsing → stability status
    ▼
Interpretation (pigean/interpretation.py)
    │ top genes, gene sets, gene-pathway links
    ▼
Report Generation (pigean/report.py)
    │ report.txt + report.html
    ▼
Provenance Manifest (pigean/manifest.py)
    │ run_manifest.json
    ▼
Output Directory (complete)
```

## Evidence Adapter Pattern

The wrapper supports multiple evidence input types through an adapter pattern:

```
BaseEvidenceAdapter (pigean/adapters/base.py)
    ├── PositiveControlsAdapter — gene lists
    ├── BayesFactorAdapter — pre-computed Bayes factors
    ├── ZscoreAdapter — Z-scores per gene
    ├── PercentileAdapter — percentile ranks per gene
    ├── ExomeAdapter — exome-based evidence
    └── GWASAdapter — GWAS summary statistics
```

Each adapter implements:
- `validate(input_path)` → file structure validation
- `normalize(input_path, output_dir)` → format normalization
- `run_qc(normalized_path, ...)` → evidence-specific QC
- `build_engine_args(normalized_path)` → engine CLI arguments

## Convergence Assessment

The wrapper assesses the PIGEAN engine's built-in stability criterion:

- **PIGEAN stability criterion**: max-fractional-SEM < 0.01 (threshold from
  `engine/priors.py --max-frac-sem`)
- **Assessment source**: `TLOG_BASED` (log file parsed) or `SUMMARY_BASED`
  (p.out only)
- **Formal MCMC convergence**: always `NOT_FORMALLY_ASSESSED` — the engine does
  not retain independent chain traces for R-hat/ESS diagnostics

See `docs/convergence.md` for full details.

## Non-Determinism

The PIGEAN engine (`engine/priors.py`) never calls `np.random.seed()`. All
stochastic outputs (gene posteriors, gene-set effects, model parameters) vary
across runs. The golden regression test uses empirically derived tolerances
from 6 independent runs to accommodate this.

## Phase History

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Golden baseline establishment | LOCKED |
| 1 | CLI wrapper + provenance | LOCKED |
| 2 | Custom gene set support | LOCKED |
| 3 | Multi-evidence adapters | LOCKED |
| 4 | GWAS summary statistics | LOCKED |
| 5 | Interpretation + stability | LOCKED |
| 5.5 | Repository reorganization | LOCKED |
