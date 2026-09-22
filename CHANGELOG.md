# Changelog

## [5.0.0] — Phase 5.5: Repository Reorganization (2026-09-21)

### Reorganized
- Moved `priors.py` → `engine/priors.py` (frozen engine clearly separated)
- Renamed `wdl/` → `workflows/`
- Renamed `ex/` → `examples/`
- Moved `gclouddock.sh` → `scripts/gclouddock.sh`
- Moved `raw/pigean.sh` → `scripts/legacy/pigean.sh`
- Moved `outputs/` → `results/development/`
- Moved `runs/` → `results/runs/`
- Collected stray root-level `.out` files → `results/development/scratch_outputs/`
- Moved `plan.README.md`, `execution.md`, `results.README.md` → `docs/phases/`
- Updated Docker tag in `scripts/gclouddock.sh` from 1.0.0 to 5.0.0

### Added
- `docs/architecture.md` — system architecture documentation
- `docs/repository_structure.md` — directory layout guide
- `docs/path_dependency_audit.md` — path dependency analysis
- `docs/repository_inventory_before.md` — pre-reorganization snapshot
- `CHANGELOG.md` — this file

### Updated
- `pigean/engine.py` — engine path updated to `engine/priors.py`
- `pigean/references.py` — Docker sentinel updated to `/app/engine/priors.py`
- `pigean/manifest.py` — priors.py path updated for SHA-256 computation
- `tests/golden/run_golden_test.sh` — engine path in fresh-run command
- `Dockerfile` — `ADD engine/priors.py` directive
- `README.md` — updated paths and added repository structure section
- `.gitignore` — added `__pycache__/`, `.pytest_cache/`, `results/` patterns

### NOT Changed
- `engine/priors.py` — byte-for-byte identical (SHA-256 verified)
- `data/*` — all reference files unchanged (checksums verified)
- Scientific behavior, statistical calculations, golden tolerances
- Test assertions or fixture data

## [5.0.0] — Phase 5: Interpretation + Stability Assessment (2026-09-21)

- PIGEAN stability assessment with max-fractional-SEM criterion
- Terminology correction: PIGEAN_STABILITY_CRITERION_MET replaces CONVERGED
- Formal MCMC convergence tracked separately (NOT_FORMALLY_ASSESSED)
- Self-contained HTML report with stability section
- Extended text report with top genes, gene sets, metric guide
- Convergence documentation (`docs/convergence.md`)

## [4.0.0] — Phase 4: GWAS Summary Statistics (2026-09-15)

- GWASAdapter for summary statistics input
- Automatic chromosome/position/effect/p-value column detection
- GWAS-specific validation and QC

## [3.0.0] — Phase 3: Multi-Evidence Adapters (2026-09-14)

- Evidence adapter architecture (BaseEvidenceAdapter)
- BayesFactorAdapter, ZscoreAdapter, PercentileAdapter, ExomeAdapter
- Evidence input documentation (`docs/evidence_inputs.md`)

## [2.0.0] — Phase 2: Custom Gene Sets (2026-09-13)

- Custom gene set loading (GMT format)
- Gene set validation and preparation
- Multiple gene set profiles (default, mouse-only, msigdb-only, custom)

## [1.0.0] — Phase 1: CLI Wrapper + Provenance (2026-09-12)

- `run_pigean.py` CLI entry point
- Input validation, normalization, gene QC
- Configuration resolution with presets
- Engine command building matching Phase 0 reference
- Run provenance manifest (`run_manifest.json`)
- Basic text report
- Golden regression test (56 checkpoints)

## [0.0.0] — Phase 0: Golden Baseline (2026-09-11)

- 6 independent engine runs for stochastic tolerance estimation
- Frozen golden reference outputs
- Golden regression test infrastructure
