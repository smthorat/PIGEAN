# Path Dependency Audit

Conducted: 2026-09-21, prior to Phase 5.5 reorganization.

## Central Path Resolution

All runtime paths flow through `pigean/references.py`:

```python
def resolve_base_dir():
    if os.path.isfile("/app/priors.py"):   # Docker sentinel
        return "/app"
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
```

`base_dir` is then used by:
- `pigean/engine.py:48` → `os.path.join(base_dir, "priors.py")`
- `pigean/references.py` → `os.path.join(base_dir, "data", ...)` for all reference files
- `pigean/manifest.py:106` → `os.path.join(base_dir, "priors.py")` for SHA-256

## `priors.py` Path Dependencies

| File | Line | Reference | Type |
|------|------|-----------|------|
| `pigean/engine.py` | 48 | `os.path.join(base_dir, "priors.py")` | Runtime path |
| `pigean/references.py` | 93 | `os.path.isfile("/app/priors.py")` | Docker sentinel |
| `pigean/manifest.py` | 106 | `os.path.join(base_dir, "priors.py")` | SHA-256 computation |
| `tests/golden/run_golden_test.sh` | 103 | `"$BASEDIR/priors.py"` | Engine invocation |
| `Dockerfile` | 12 | `ADD priors.py priors.py` | Docker image build |
| `tests/test_engine_command.py` | 38 | `endswith("priors.py")` | Test assertion (suffix-based — safe) |

## `data/` Path Dependencies

| File | Lines | Count | Reference Pattern |
|------|-------|-------|-------------------|
| `pigean/references.py` | 24-50 | ~10 | `os.path.join("data", "filename")` |
| `Dockerfile` | 5-10 | 6 | `ADD data/filename data/` |
| `tests/golden/run_golden_test.sh` | 104-106+ | 6 | `"$BASEDIR/data/filename"` |
| `tests/test_engine_command.py` | 18-26 | 6 | `os.path.join(REPO_ROOT, "data", ...)` |

## `sys.path` Dependencies

All test files and `run_pigean.py` insert repo root into `sys.path`:

| File | Pattern |
|------|---------|
| `run_pigean.py:10` | `sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))` |
| All `tests/test_*.py:3` | `sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))` |

These resolve to repo root because `run_pigean.py` is at root and tests are in `tests/`. Moving either would break imports.

## Docker Path Dependencies

```dockerfile
WORKDIR /app
ADD data/...        data/          # data/ relative to /app
ADD priors.py       priors.py      # engine at /app/priors.py
ADD run_pigean.py   run_pigean.py  # entry point at /app/run_pigean.py
ADD pigean/         pigean/        # package at /app/pigean/
ADD docs/           docs/          # docs at /app/docs/
ENTRYPOINT ["python3", "-u", "/app/run_pigean.py"]
```

## WDL Path Dependencies

`wdl/rock_pigean.wdl` command line: `python3 -u /app/run_pigean.py \` — references the Docker container path, not the repo path. Renaming the `wdl/` directory has no effect on WDL execution.

## Risk Assessment for Moves

| Move | Risk | Path Changes Required |
|------|------|-----------------------|
| `priors.py` → `engine/priors.py` | LOW | 3 Python + 1 bash + 1 Dockerfile |
| `data/` → `resources/` | MODERATE | ~28 path references across 4 files |
| `wdl/` → `workflows/` | NONE | Directory rename only; no code references |
| `ex/` → `examples/` | NONE | README text only |
| `outputs/` → `results/development/` | NONE | No code references |
| `runs/` → `results/runs/` | NONE | No code references |
| `raw/` → `scripts/legacy/` | NONE | No code references |
| `pigean/` → anywhere | HIGH | Breaks all imports |
| `run_pigean.py` → anywhere | HIGH | Docker, WDL, sys.path |
| `tests/` → anywhere | MODERATE | sys.path, pytest discovery |
