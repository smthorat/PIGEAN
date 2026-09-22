# Implementation Plan: PIGEAN Wrapper Layer (v1 — Phases 0–2)

## Context

The rock-pigean repo (`broadinstitute/rock-pigean`) contains a 16,567-line PIGEAN engine (`priors.py`) with 309 `optparse` flags, but the WDL workflow exposes only 2 inputs and hardcodes everything else. Scientists can't customize gene sets, validate their inputs, or understand outputs without manual `awk`. The plan.README.md describes a phased generalization — this plan implements Phases 0–2 (the v1 release: golden test + wrapper + custom gene sets) inside the same repo, without touching `priors.py`.

**Decision made:** All new code lives in the same repo alongside existing files.

---

## Target File Structure After v1

```
rock-pigean/
│
├── run_pigean.py                  ← NEW: CLI entry point
├── pigean/                        ← NEW: Python package
│   ├── __init__.py                ←   version string
│   ├── config.py                  ←   defaults, presets, config merge (single source of truth)
│   ├── validation.py              ←   input gene QC, background QC, combined status
│   ├── engine.py                  ←   build + run priors.py command
│   ├── manifest.py                ←   run_manifest.json generation (always written, even on FAIL)
│   ├── report.py                  ←   basic text report (verified metrics only)
│   ├── references.py              ←   reference/annotation path resolution
│   └── adapters/
│       ├── __init__.py            ←   base adapter class + registry
│       ├── positive_controls.py   ←   gene-list adapter (v1)
│       └── gene_sets.py           ←   gene set loading, GMT conversion (v1 Phase 2)
│
├── priors.py                      ← UNCHANGED (engine-freeze)
├── Dockerfile                     ← UPDATED: adds run_pigean.py + pigean/
├── wdl/rock_pigean.wdl            ← UPDATED: parameterized, calls run_pigean.py
├── data/                          ← UNCHANGED (reference files stay here for v1)
├── ex/gene_list                   ← UNCHANGED
├── raw/pigean.sh                  ← UNCHANGED
│
├── tests/                         ← NEW: frozen test fixtures + unit tests (checked in)
│   ├── golden/
│   │   ├── input/
│   │   │   └── positive_controls.txt
│   │   ├── expected/
│   │   │   ├── gs.golden.out
│   │   │   ├── gss.golden.out
│   │   │   ├── ggss.golden.out
│   │   │   └── p.golden.out
│   │   ├── checkpoints.json
│   │   └── run_golden_test.sh
│   ├── test_validation.py
│   ├── test_config.py             ← tests config precedence (defaults < JSON < CLI)
│   ├── test_gene_sets.py          ← tests GMT conversion, replace vs supplement
│   └── test_engine_command.py     ← proves wrapper generates correct effective engine settings
│
├── outputs/                       ← NEW: development validation artifacts (gitignored)
│   ├── phase0_golden/
│   ├── phase1_wrapper/
│   ├── phase2_gene_sets/
│   ├── phase3_evidence/
│   ├── phase4_gwas/
│   ├── phase5_report/
│   ├── phase6_advanced/
│   └── phase7_production/
│
├── runs/                          ← NEW: normal researcher analyses (gitignored)
│
└── README.md                      ← UPDATED: usage docs
```

### Key Structural Decisions

- **`tests/golden/`** = frozen reference fixtures (checked into git). Small, static, the regression baseline.
- **`outputs/phaseX/`** = generated development artifacts (gitignored). Your audit trail — what you ran, what changed, what you compared. Never checked in.
- **`runs/`** = normal researcher analyses once the pipeline is mature (gitignored).
- **`data/`** stays as-is for v1. No reorganization into `references/hg19/` yet — that's Phase 4.
- **`run_pigean.py`** is a thin CLI entry point (~60 lines). All logic lives in `pigean/`.
- **No `setup.py`/`pyproject.toml` for v1** — the code runs from the repo directory inside Docker.
- **No new dependencies** — use JSON (stdlib) instead of YAML for config/manifest files. Everything the wrapper needs is Python stdlib + what `priors.py` already requires (scipy/numpy).
- **`config.py` is the single source of truth for defaults** — argparse uses `default=None` so that only explicitly supplied CLI values override the config.
- **File paths come from CLI/WDL; settings come from config.** The JSON config file (`--config`) contains parameters and settings (gene_sets, genome_build, preset, custom_gene_set_action, etc.). Actual input file paths (primary input, custom gene-set files, background) always come from the CLI or WDL — never from the JSON config. This keeps the design simple and avoids ambiguity about which layer owns file resolution.

### Universal Run Output Pattern

Every run — development or production — writes this standard layout into its output directory. The **run directory itself identifies the experiment**, so output filenames are simple (`gs.out`, not `gs.<name>.out`):

```
<run_dir>/
├── input/
│   └── <original input file(s)>          ← what the user provided
├── normalized/
│   └── <actual input sent to engine>     ← what the wrapper transformed
├── priors_command.txt                     ← exact engine command used
├── resolved_config.json                   ← full merged config (defaults + overrides)
├── input_gene_qc.tsv                      ← per-gene validation status
├── run_manifest.json                      ← provenance — ALWAYS written, even on FAIL
├── pigean_run.log                         ← captured engine stdout/stderr
├── gs.out                                 ← gene stats
├── gss.out                                ← gene set stats
├── ggss.out                               ← gene × gene set stats
├── p.out                                  ← model parameters
└── report.txt                             ← basic text report
```

This answers three questions for every run: **What went in? What was changed? What came out?**

> **Failed runs still leave an audit trail.** If validation fails, the run directory still contains `run_manifest.json` (with `execution_status: "NOT_RUN"`), `resolved_config.json`, and `input_gene_qc.tsv` (if QC was reached). The manifest explains *why* nothing else exists.

> **Output directory must be empty.** If the output directory already contains results from a previous run, the wrapper fails with a clear error rather than risk mixing old and new outputs. Pass `--overwrite` to force — this **deletes all old contents** and recreates a clean directory, so no stale files can survive from the previous run.

### Correct Pipeline Execution Order

```
1. Create CLEAN run directory      (fail if non-empty without --overwrite;
                                    --overwrite deletes old contents first)
         ↓
2. Resolve configuration           (DEFAULTS → preset → config JSON → explicit CLI)
         ↓
3. Save resolved_config.json       (present even for the earliest failures)
         ↓
4. Validate primary input file     (file exists, readable, non-empty, looks like text)
         ↓  FAIL → write manifest (execution_status: NOT_RUN) + exit
         ↓
5. Preserve original input         (copy to input/)
         ↓
6. Normalize primary input         (strip whitespace, skip blanks → normalized/)
         ↓
7. If background requested:
      validate background file
         ↓  FAIL → write manifest (execution_status: NOT_RUN) + exit immediately
         ↓
      normalize background file    (→ input/background_original.txt + normalized/background.txt)
         ↓
8. Validate + prepare custom gene sets (if provided)
      validate files exist, format/action specified
      convert GMT → engine format into normalized/
         ↓  FAIL → write manifest (execution_status: NOT_RUN) + exit
         ↓         (short-circuit BEFORE path resolution)
         ↓
9. Resolve final gene-set paths + reference paths
      uses NORMALIZED custom paths (not originals)
         ↓
10. Gene / reference QC            (check NORMALIZED genes against gene map, loc files, gene sets)
         ↓
11. Background QC                  (if background provided: consistency with input genes)
         ↓
12. Combined validation            (merge all stages → PASS / WARN / FAIL)
         ↓  FAIL → write QC + manifest (execution_status: NOT_RUN) + exit
         ↓
13. Build priors.py command
         ↓
14. Save priors_command.txt
         ↓
15. Run PIGEAN                     (subprocess, capture log)
         ↓
16. gs / gss / ggss / p
         ↓
17. Final manifest + verified report
```

> **Important:** Gene QC runs on the **normalized** input, not the raw file. This ensures `input_gene_qc.tsv` describes the actual genes sent to `priors.py`.

### Phase Names Are Development Tracking Only

The production code (`run_pigean.py`, `pigean/`) knows only scientific concepts: analysis type, input type, gene-set profile, background, genome build, preset, output directory. It never references "phase0", "phase1", etc. Phase folders exist solely for organizing your development experiments.

---

## Phase 0 — Golden Test (~0.5 day)

### What
Capture the current known-good run as a regression baseline. Every future change proves it didn't break science.

### Output Structure

```
outputs/phase0_golden/
├── run_01/
│   ├── input/
│   │   └── positive_controls.txt
│   ├── priors_command.txt
│   ├── gs.out
│   ├── gss.out
│   ├── ggss.out
│   ├── p.out
│   └── pigean_run.log
├── run_02/
├── run_03/
├── run_04/
├── run_05/
│
├── baseline_summary/
│   ├── variability_summary.tsv        ← per-gene prior mean/stddev across runs
│   ├── top_gene_overlap.tsv           ← top-N gene agreement across runs
│   ├── top_pathway_overlap.tsv        ← top-N pathway agreement across runs
│   ├── rank_correlations.tsv          ← pairwise Spearman ρ between runs
│   └── final_tolerances.json          ← calibrated acceptance thresholds
│
└── summary.txt                         ← purpose, result, open issues
```

### Steps

1. **Create `tests/golden/input/positive_controls.txt`** — copy of `ex/gene_list`

2. **Create `tests/golden/expected/`** — copy current `gs.out`, `gss.out`, `ggss.out`, `p.out` as frozen golden references

3. **Create `tests/golden/checkpoints.json`** — key values to validate. **All tolerances are provisional until Phase 0 repeated runs determine which quantities are deterministic.**
   ```json
   {
     "_note": "Tolerances are provisional. Phase 0 repeated runs determine which quantities are deterministic (exact checks) vs stochastic (empirical tolerances).",
     "structural": {
       "_note": "Assumed deterministic — Phase 0 will confirm or add tolerances.",
       "gs_rows": {"expected": 42216, "tolerance": "TBD"},
       "gss_rows": {"expected": 41627, "tolerance": "TBD"},
       "p_rows": {"expected": 36, "tolerance": "TBD"}
     },
     "model_parameters": {
       "num_genes_read": 42216,
       "num_gene_sets_read": 2088,
       "num_chains": 10,
       "num_gibbs_iter": {
         "_note": "p.out reports 499. Confirm whether this is the final zero-based iteration index (500 iterations performed) or the actual count. Store both once clarified.",
         "reported_value": 499
       }
     },
     "gene_priors": {
       "LEP":    {"expected": 1.89,    "tolerance": 0.1},
       "ADIPOQ": {"expected": 1.36,    "tolerance": 0.1},
       "LEPR":   {"expected": 0.988,   "tolerance": 0.1},
       "PLINK1": {"expected": -0.0104, "tolerance": 0.05}
     },
     "rank_stability": {
       "top5_non_input_must_include": ["LEPR", "PPARG", "APOE"]
     }
   }
   ```

4. **Create `tests/golden/run_golden_test.sh`** — bash script that validates using **tolerance-based comparison** (not byte-for-byte diff):
   - For each quantity, Phase 0 repeated runs determine the check type:
     - **Deterministic** (identical across all 5 runs) → exact match
     - **Stochastic** (varies across runs) → empirical tolerance from observed range
   - Structural counts → likely deterministic (verify)
   - Key gene priors → tolerance from observed variability
   - Top-N gene overlap → threshold
   - Gene ranking correlation → Spearman ρ threshold
   - Top-N pathway overlap → threshold
   - Reports PASS/FAIL

5. **Run the engine 5 times** into `outputs/phase0_golden/run_01` through `run_05` to measure stochastic variability (the engine seeds `random.seed(0)` but may not seed `np.random`).

6. **Generate `baseline_summary/`** — compute variability metrics across the 5 runs, calibrate final tolerances, and update `checkpoints.json` with empirical thresholds.

> **Environment note:** Record the exact Docker image ID/digest used for the Phase 0 baseline runs. The Dockerfile uses `ubuntu:plucky` plus `apt upgrade`, so the environment can change between builds. Container pinning is Phase 7, but knowing the baseline environment is essential for interpreting any future variability.

---

## Phase 1 — Wrapper Layer + Parameterized WDL (~4–5 days)

### Output Structure

```
outputs/phase1_wrapper/
├── original_command/
│   ├── priors_command.txt                ← the exact WDL command (our specification)
│   ├── gs.out
│   ├── gss.out
│   ├── ggss.out
│   └── p.out
│
├── wrapper_command/
│   ├── input/
│   │   └── positive_controls.txt
│   ├── normalized/
│   │   └── positive_controls.txt
│   ├── priors_command.txt                ← wrapper-generated command
│   ├── resolved_config.json
│   ├── input_gene_qc.tsv
│   ├── run_manifest.json
│   ├── pigean_run.log
│   ├── gs.out
│   ├── gss.out
│   ├── ggss.out
│   ├── p.out
│   └── report.txt
│
├── comparison/
│   ├── command_diff.txt                  ← diff of the two priors_command.txt files (informational)
│   ├── parameter_comparison.txt          ← effective parameters side-by-side (the real test)
│   ├── gene_score_comparison.tsv
│   ├── rank_correlation.txt
│   ├── pathway_overlap.tsv
│   └── compatibility_result.txt          ← PASS / FAIL
│
└── summary.txt
```

The central question Phase 1 answers: **Did the wrapper change the science?** The `comparison/` folder makes that directly inspectable.

> **Note on command comparison:** The wrapper may legitimately produce different paths (e.g., `normalized/positive_controls.txt` vs `ex/gene_list`). The PASS criterion is that the **effective engine parameters and reference files are identical**, not that the raw command strings match character-for-character. `command_diff.txt` is saved for informational debugging; `parameter_comparison.txt` is the actual test.

### Implementation Order

Build in this sequence — each step is testable against the golden baseline:

#### Step 1: `pigean/config.py` — Defaults and Config Merge

**Responsibility:** Define all default values, presets, and the config-merge logic. `config.py` is the **single source of truth** for defaults — argparse does not define its own.

```python
# pigean/config.py

DEFAULTS = {
    "analysis": "positive-controls",
    "preset": "standard",
    "genome_build": "hg19",
    "gene_sets": "default",
    "max_num_gene_sets": 5000,
    "gene_filter_value": 1,
    "gene_set_filter_value": 0.01,
    "debug_level": 3,
    "num_chains": 10,       # implicit in engine default
}

# Paths relative to /app inside the Docker container
REFERENCE_PATHS = {
    "hg19": {
        "gene_loc": "data/NCBI37.3.plink.gene.loc",
        "tss_loc": "data/refGene_hg19_TSS.subset.loc",
        "exons_loc": "data/NCBI37.3.plink.gene.exons.loc",
    }
}

ANNOTATION_PATHS = {
    "default": [
        "data/gene_set_list_mouse_2024.txt",
        "data/gene_set_list_msigdb_nohp.txt",
    ],
    "mouse-only": ["data/gene_set_list_mouse_2024.txt"],
    "msigdb-only": ["data/gene_set_list_msigdb_nohp.txt"],
}

GENE_MAP_PATH = "data/portal_gencode.gene.map"

def resolve_config(args, config_file=None):
    """Merge in strict precedence order:
    
        DEFAULTS  →  preset overrides  →  config JSON file  →  explicit CLI args
    
    argparse defaults are all None, so only values the user actually typed
    on the command line override the config file. This prevents argparse
    defaults from silently clobbering JSON settings.
    
    Example: if config.json says gene_sets=msigdb-only and the user didn't
    pass --gene-sets on the CLI, the result is msigdb-only (not "default").
    """
    config = dict(DEFAULTS)
    # 1. Apply preset overrides (if any future presets change defaults)
    # 2. Load and merge config_file (JSON) if provided
    if config_file:
        with open(config_file) as f:
            overrides = json.load(f)
        config.update(overrides)
    # 3. Apply only CLI args that were explicitly supplied (not None)
    for key in [
        "analysis",
        "gene_sets",
        "genome_build",
        "preset",
        "custom_gene_set_action",
        "custom_gene_set_format",
    ]:
        cli_val = getattr(args, key.replace("-", "_"), None)
        if cli_val is not None:
            config[key] = cli_val
    # Record background and custom files in config (for resolved_config.json)
    if getattr(args, "background", None):
        config["background_provided"] = True
    if getattr(args, "custom_gene_set_files", None):
        config["custom_gene_set_files"] = args.custom_gene_set_files
    return config

def write_resolved_config(config, output_dir):
    """Write resolved_config.json to the output directory.
    
    Includes all scientific choices so you can inspect one file and
    understand the entire run configuration:
    
        {
          "analysis": "positive-controls",
          "preset": "standard",
          "genome_build": "hg19",
          "gene_sets": "msigdb-only",
          "custom_gene_set_format": "gmt",
          "custom_gene_set_action": "supplement",
          "background_provided": true,
          "max_num_gene_sets": 5000,
          ...
        }
    """
    ...
```

#### Step 2: `pigean/references.py` — Path Resolution

**Responsibility:** Resolve symbolic names ("hg19", "default") to actual file paths. Handles both in-Docker (`/app/data/...`) and local-dev paths. Also assembles the final `--X-in` list from base annotations + prepared custom paths.

> **Single responsibility note:** `gene_sets.py` validates and converts custom gene-set files (returning `GeneSetPreparationResult` with normalized paths). `references.py::resolve_gene_set_paths()` assembles the final list of `--X-in` paths by combining base annotations with those prepared paths. There is exactly one function that creates the final engine gene-set list.

```python
# pigean/references.py

def resolve_base_dir():
    """Return /app if running in Docker, else the repo root."""
    if os.path.exists("/app/priors.py"):
        return "/app"
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def resolve_reference_paths(genome_build, base_dir):
    """Return dict of {gene_loc, tss_loc, exons_loc} absolute paths."""
    ...

def resolve_gene_set_paths(gene_set_profile, custom_normalized_paths,
                           custom_action, base_dir):
    """Return list of gene set file paths ready for the engine.
    
    custom_normalized_paths are ALREADY validated and converted (by
    prepare_custom_gene_sets). This function only assembles the final
    list based on the profile and action — it does no file I/O.
    
    Handles replace vs supplement against the base profile."""
    ...

def resolve_gene_map_path(base_dir):
    """Return absolute path to gene map file."""
    ...
```

#### Step 3: `pigean/adapters/positive_controls.py` — First Adapter

**Responsibility:** Basic file validation, normalization, and engine-arg construction for gene-list input.

```python
# pigean/adapters/positive_controls.py
from pigean.adapters import EvidenceAdapter, FileValidationResult, NormalizationResult

class PositiveControlsAdapter(EvidenceAdapter):
    def validate_file(self, input_path):
        """Basic file checks: exists, readable, non-empty, looks like text with one gene per line.
        This runs BEFORE normalization and BEFORE copying to input/."""
        # Return FileValidationResult(status=PASS|FAIL, issues=[...])
    
    def normalize(self, input_path, output_dir, label="positive_controls"):
        """Strip whitespace, skip blank lines, deduplicate.
        Copies original to input/<label>.txt, writes cleaned to normalized/<label>.txt.
        This runs AFTER validate_file but BEFORE gene QC.
        
        label="positive_controls" → input/positive_controls.txt, normalized/positive_controls.txt
        label="background"        → input/background_original.txt, normalized/background.txt
        """
        # Return NormalizationResult(normalized_path=..., input_count=..., ...)
    
    def build_engine_args(self, norm_result, config):
        """Return ['--positive-controls-in', '<normalized_path>']
        
        Called by engine.py — adapter owns the translation from
        NormalizationResult (data) to engine flags (behavior)."""
        return ["--positive-controls-in", norm_result.normalized_path]
```

**Base adapter contract** in `pigean/adapters/__init__.py`:
```python
from dataclasses import dataclass, field
from enum import Enum

class ValidationStatus(Enum):
    PASS = "PASS"
    PASS_WITH_WARNINGS = "PASS_WITH_WARNINGS"
    FAIL = "FAIL"

@dataclass
class FileValidationResult:
    """Result of basic file-level checks (before normalization)."""
    status: ValidationStatus
    issues: list = field(default_factory=list)

@dataclass
class NormalizationResult:
    """Result of normalization (after file validation, before gene QC).
    This is data — it does not own behavior (the adapter does)."""
    normalized_path: str
    input_count: int = 0
    excluded: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

@dataclass
class GeneQCResult:
    """Result of gene/reference QC (after normalization)."""
    status: ValidationStatus
    input_count: int = 0
    recognized: int = 0
    with_coordinates: int = 0
    with_annotations: int = 0
    unresolved: list = field(default_factory=list)
    per_gene_records: list = field(default_factory=list)

@dataclass
class GeneSetPreparationResult:
    """Result of custom gene-set validation + conversion.
    Carries the NORMALIZED paths that the engine should receive."""
    status: ValidationStatus
    original_paths: list = field(default_factory=list)
    normalized_paths: list = field(default_factory=list)
    issues: list = field(default_factory=list)

class EvidenceAdapter:
    def validate_file(self, input_path): ...
    def normalize(self, input_path, output_dir, label="input"): ...
    def build_engine_args(self, norm_result, config): ...
```

#### Step 4: `pigean/validation.py` — Gene QC + Background QC + Combined Status

**Responsibility:** Load reference files, check each **normalized** input gene against them, check background consistency, produce `input_gene_qc.tsv`, and compute the **combined validation status**.

```python
# pigean/validation.py

def load_gene_map(gene_map_path):
    """Load portal_gencode.gene.map → set of known gene symbols."""
    ...

def load_gene_locations(loc_file_path):
    """Load gene.loc → dict of {gene_symbol: (chrom, start, end)}."""
    ...

def count_gene_set_memberships(gene, gene_set_files):
    """Count how many gene sets contain this gene across all gene set files."""
    ...

def validate_normalized_genes(normalized_genes, gene_map, gene_locs, gene_set_files):
    """Run gene/reference QC on the NORMALIZED gene list (not the raw input).
    Returns GeneQCResult."""
    ...

def validate_background(background_norm_result, input_norm_result, gene_map):
    """Check background gene consistency:
    - How many background genes resolve in gene map?
    - Are ALL input genes present in the background?
    - Report any input genes missing from background.
    
    Returns BackgroundQCResult with:
      total, mapped, unmapped, input_genes_missing_from_background"""
    ...

def prepare_custom_gene_sets(custom_files, custom_format, custom_action,
                             profile, output_dir):
    """Validate gene-set configuration and return NORMALIZED paths.
    
    Called UNCONDITIONALLY (even with no custom files) so the full
    legal-combinations table is enforced in one place:
    
      default/mouse-only/msigdb-only + no custom  → PASS (no-op)
      custom + no files                            → FAIL
      files + no format                            → FAIL
      files + no action                            → FAIL
      action + no files                            → FAIL
      files + format + action                      → validate, convert, PASS
    
    All values come from the RESOLVED config (not raw CLI args), so
    JSON config overrides are respected.
    
    Steps (when custom files are provided):
      1. Check custom_format is provided (required with files)
      2. Check custom_action is provided (required with files)
      3. Validate all custom files exist and are readable
      4. If format=gmt: convert each to engine format in normalized/
         If format=engine: copy to normalized/ as-is
      5. Warn if gene sets have very few entries
    
    Returns GeneSetPreparationResult:
      - original_paths:   what the user provided (for input/ and manifest)
      - normalized_paths:  what the engine should receive (converted if GMT)
      - status:           PASS / PASS_WITH_WARNINGS / FAIL
      - issues:           list of problems found
    
    Audit trail:
      input/custom_pathways.gmt            ← original
      normalized/custom_pathways.engine.txt ← converted (what PIGEAN receives)
    """
    ...

def compute_combined_status(file_result, norm_result, gene_qc_result,
                            background_file_result=None,
                            background_qc_result=None,
                            custom_gs_prep_result=None):
    """Merge all validation stages into one final status.
    
    FAIL if:
      - file validation failed
      - 0 recognized genes after normalization
      - no usable gene sets loaded
      - background was REQUESTED but its file validation failed
        (a requested-but-invalid background must never be silently ignored)
      - custom gene-set preparation failed
        (missing file, missing format/action, conversion error)
    
    PASS_WITH_WARNINGS if:
      - some genes unresolved but at least 1 usable
      - normalization removed duplicates or blank lines
      - input genes missing from declared background
    
    PASS if:
      - all checks passed, all genes resolved
    
    Returns (ValidationStatus, list of all issues from all stages)"""
    ...

def write_gene_qc_tsv(gene_qc_result, output_dir):
    """Write input_gene_qc.tsv."""
    ...
```

#### Step 5: `pigean/engine.py` — Build and Run priors.py

**Responsibility:** Construct the exact `priors.py` command line from resolved config, execute it, and save the command to `priors_command.txt`.

> **Engine flag note:** The engine defines `--gene-gene-set-stats-out` (priors.py line 228). The current WDL uses the prefix shorthand `--gene-gene-set-stats` which works via optparse prefix matching. The wrapper uses the **full flag name** for clarity.

> **Adapter/builder responsibility split:** The adapter owns `build_engine_args(norm_result, config)` → returns evidence-specific flags. The engine builder owns everything else (references, filters, output paths). `NormalizationResult` is data; the adapter is behavior.

```python
# pigean/engine.py

def build_priors_command(config, adapter, norm_result, reference_paths,
                         gene_set_paths, gene_map_path, output_dir,
                         background_norm_result=None):
    """
    Build the priors.py command list.
    
    Responsibilities:
      - adapter.build_engine_args(norm_result, config) → evidence-specific flags
      - this function → everything else (references, filters, outputs)
    
    Uses full flag names (not optparse prefix shortcuts):
        python3 -u /app/priors.py gibbs \
            --X-in <gene_set_1> --X-in <gene_set_2> \
            --gene-map-in <gene_map> \
            --max-num-gene-sets 5000 \
            --gene-stats-out <output_dir>/gs.out \
            --gene-set-stats-out <output_dir>/gss.out \
            --gene-gene-set-stats-out <output_dir>/ggss.out \
            --params-out <output_dir>/p.out \
            --debug-level 3 \
            --positive-controls-in <normalized_input> \
            --gene-loc-file <gene.loc> \
            --gene-loc-file-huge <tss.loc> \
            --exons-loc-file-huge <exons.loc> \
            --gene-filter-value 1 --gene-set-filter-value 0.01
    """
    cmd = ["python3", "-u", os.path.join(base_dir, "priors.py"), "gibbs"]
    for gs_path in gene_set_paths:
        cmd.extend(["--X-in", gs_path])
    cmd.extend(["--gene-map-in", gene_map_path])
    cmd.extend(["--max-num-gene-sets", str(config["max_num_gene_sets"])])
    cmd.extend(["--gene-stats-out", os.path.join(output_dir, "gs.out")])
    cmd.extend(["--gene-set-stats-out", os.path.join(output_dir, "gss.out")])
    cmd.extend(["--gene-gene-set-stats-out", os.path.join(output_dir, "ggss.out")])
    cmd.extend(["--params-out", os.path.join(output_dir, "p.out")])
    cmd.extend(["--debug-level", str(config["debug_level"])])
    # Evidence-specific flags from the adapter (behavior on data)
    cmd.extend(adapter.build_engine_args(norm_result, config))
    cmd.extend(["--gene-loc-file", reference_paths["gene_loc"]])
    cmd.extend(["--gene-loc-file-huge", reference_paths["tss_loc"]])
    cmd.extend(["--exons-loc-file-huge", reference_paths["exons_loc"]])
    cmd.extend(["--gene-filter-value", str(config["gene_filter_value"])])
    cmd.extend(["--gene-set-filter-value", str(config["gene_set_filter_value"])])
    # Background uses its own normalized path (not the raw config path)
    if background_norm_result:
        cmd.extend(["--positive-controls-all-in", background_norm_result.normalized_path])
    return cmd

def save_command(cmd, output_dir):
    """Write the exact priors.py command to priors_command.txt for debugging and comparison."""
    with open(os.path.join(output_dir, "priors_command.txt"), "w") as f:
        f.write(" \\\n    ".join(cmd) + "\n")

def run_engine(cmd, output_dir, log_file="pigean_run.log"):
    """Execute priors.py, capture stdout/stderr to log file, return exit code."""
    save_command(cmd, output_dir)
    ...
```

#### Step 6: `pigean/manifest.py` — Run Manifest

The manifest is **always written**, including for failed validation runs (`execution_status: "NOT_RUN"`). This ensures every attempted run has provenance.

```python
# pigean/manifest.py

def generate_manifest(config, combined_status, all_issues, gene_qc,
                      engine_exit_code, output_dir, start_time, end_time,
                      input_path=None, normalized_path=None,
                      gene_map_path=None, gene_set_paths=None,
                      reference_paths=None,
                      background_qc=None, background_path=None,
                      background_norm_path=None):
    """Write run_manifest.json with full provenance.
    
    ALWAYS called — even when validation fails (engine_exit_code=None).
    Uses the COMBINED validation status (file + normalization + gene QC + background QC).
    Records ALL scientific choices so one file explains the entire run configuration."""
    manifest = {
        "wrapper_version": pigean.__version__,
        "engine": {
            "name": "priors.py",
            "sha256": compute_file_sha256(priors_py_path) if priors_py_path else None,
        },
        "dependencies": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "analysis_type": config["analysis"],
        "preset": config["preset"],
        "genome_build": config["genome_build"],
        "gene_sets": config["gene_sets"],
        
        # Custom gene-set configuration
        "custom_gene_set_format": config.get("custom_gene_set_format"),
        "custom_gene_set_action": config.get("custom_gene_set_action"),
        "custom_gene_set_files": config.get("custom_gene_set_files"),
        
        # Input checksums
        "input_checksum": compute_file_sha256(input_path) if input_path else None,
        "normalized_input_checksum": compute_file_sha256(normalized_path) if normalized_path else None,
        
        # Reference checksums
        "gene_map_checksum": compute_file_sha256(gene_map_path) if gene_map_path else None,
        "gene_set_checksums": {os.path.basename(p): compute_file_sha256(p)
                               for p in (gene_set_paths or [])},
        "reference_checksums": {k: compute_file_sha256(v)
                                for k, v in (reference_paths or {}).items()},
        
        # Gene QC
        "input_genes": gene_qc.input_count if gene_qc else None,
        "recognized_genes": gene_qc.recognized if gene_qc else None,
        "unresolved_genes": gene_qc.unresolved if gene_qc else None,
        
        # Background provenance
        "background_provided": background_path is not None,
        "background_checksum": compute_file_sha256(background_path) if background_path else None,
        "background_normalized_checksum": compute_file_sha256(background_norm_path) if background_norm_path else None,
        "background_qc": {
            "total": background_qc.total,
            "mapped": background_qc.mapped,
            "unmapped": background_qc.unmapped,
            "input_genes_missing_from_background": background_qc.input_genes_missing,
        } if background_qc else None,
        
        # Validation
        "validation_status": combined_status.value,
        "validation_issues": all_issues,
        
        # Execution
        "execution_status": (
            "NOT_RUN" if engine_exit_code is None
            else "COMPLETED" if engine_exit_code == 0
            else "FAILED"
        ),
        "exit_code": engine_exit_code,
        "convergence_status": None,  # Phase 5
        
        # Timing
        "start_time": start_time,
        "end_time": end_time,
        "duration_seconds": (end_time - start_time) if (start_time and end_time) else None,
        
        "log_file": "pigean_run.log" if engine_exit_code is not None else None,
    }
    with open(os.path.join(output_dir, "run_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
```

#### Step 7: `pigean/report.py` — Basic Text Report

> **Interpretation contract prerequisite:** Every PIGEAN metric displayed in the report must first be entered into an interpretation contract with: source column, scale, transformation, and biological meaning. If a metric (e.g., `combined_D`) is not yet fully verified against the engine version, it is **not shown** in v1. The v1 report displays only verified metrics.

```python
# pigean/report.py

# INTERPRETATION CONTRACT
# -----------------------
# To be populated ONLY after tracing each field against the exact priors.py version.
# Add entries one by one as they are verified. Unverified metric = don't show it.
#
# Verified for v1:
#   (to be filled during Phase 1 implementation by tracing priors.py output)
#
# Not yet verified (deferred to Phase 5):
#   combined_D — labeled "posterior probability" but exact derivation not yet traced
#   combined   — labeled "combined score" but prior + log_bf arithmetic not fully verified
#
# Rule: no field appears in the report without an entry in this contract.

def generate_text_report(config, gene_qc, model_params, gs_path, gss_path,
                         ggss_path, output_dir):
    """
    Generate a basic text report with:
    - Input QC summary
    - Model status (completed, not assessed for convergence)
    - Ranked input genes (with verified metrics only)
    - Top 10 non-input candidate genes (with verified metrics only)
    - Top 10 enriched pathways (with verified metrics only)
    - Score definitions for every displayed field
    
    Does NOT display any metric until it has been traced against priors.py
    and entered into the interpretation contract above.
    """
    ...
```

#### Step 8: `run_pigean.py` — CLI Entry Point

```python
#!/usr/bin/env python3
"""PIGEAN Pipeline Wrapper — user-facing interface to priors.py."""

import argparse
import sys
import os
import shutil
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pigean.config import resolve_config, write_resolved_config
from pigean.references import (resolve_base_dir, resolve_reference_paths,
                                resolve_gene_set_paths, resolve_gene_map_path)
from pigean.validation import (validate_normalized_genes, validate_background,
                                prepare_custom_gene_sets,
                                compute_combined_status, write_gene_qc_tsv,
                                load_gene_map, load_gene_locations)
from pigean.adapters import ValidationStatus
from pigean.adapters.positive_controls import PositiveControlsAdapter
from pigean.engine import build_priors_command, run_engine
from pigean.manifest import generate_manifest
from pigean.report import generate_text_report

def main():
    parser = argparse.ArgumentParser(description="PIGEAN Pipeline")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--overwrite", action="store_true",
                        help="Delete old contents if output directory is non-empty")
    # All defaults are None — config.py is the single source of truth
    parser.add_argument("--analysis", default=None,
                        choices=["positive-controls"])
    parser.add_argument("--gene-sets", default=None,
                        choices=["default", "mouse-only", "msigdb-only", "custom"])
    parser.add_argument("--genome-build", default=None, choices=["hg19"])
    parser.add_argument("--preset", default=None, choices=["standard"])
    # Phase 2:
    parser.add_argument("--custom-gene-set-files", nargs="+")
    parser.add_argument("--custom-gene-set-format", choices=["engine", "gmt"],
                        help="Required when --custom-gene-set-files is used")
    parser.add_argument("--custom-gene-set-action", default=None,
                        choices=["replace", "supplement"],
                        help="Whether custom gene sets replace or supplement the base profile")
    parser.add_argument("--background")
    # Advanced:
    parser.add_argument("--config",
                        help="Path to advanced config JSON file (merged with defaults)")
    args = parser.parse_args()
    
    start_time = time.time()

    # ── Step 1: Create CLEAN run directory ──────────────────────────
    if os.path.exists(args.output) and os.listdir(args.output):
        if not args.overwrite:
            print(f"ERROR: Output directory is not empty: {args.output}")
            print("Use --overwrite to force, or choose a new directory.")
            sys.exit(1)
        # --overwrite: delete old contents so no stale files survive
        shutil.rmtree(args.output)
    os.makedirs(args.output, exist_ok=True)

    # ── Step 2–3: Resolve config and save immediately ───────────────
    # (resolved_config.json is present even for the earliest failures)
    config = resolve_config(args, config_file=args.config)
    write_resolved_config(config, args.output)

    # ── Step 4: Validate primary input file ─────────────────────────
    adapter = PositiveControlsAdapter()
    file_result = adapter.validate_file(args.input)
    if file_result.status == ValidationStatus.FAIL:
        generate_manifest(config, ValidationStatus.FAIL,
                          file_result.issues, None, None, args.output,
                          start_time, time.time())
        print_failure(file_result)
        sys.exit(1)
    
    # ── Step 5–6: Preserve original + normalize ─────────────────────
    norm_result = adapter.normalize(args.input, args.output,
                                     label="positive_controls")
    
    # ── Step 7: Background (if requested) ───────────────────────────
    # A requested-but-invalid background is a hard FAIL — exit immediately,
    # consistent with custom gene-set failure behavior.
    background_norm = None
    background_file_result = None
    background_qc = None
    if args.background:
        bg_adapter = PositiveControlsAdapter()
        background_file_result = bg_adapter.validate_file(args.background)
        if background_file_result.status == ValidationStatus.FAIL:
            # Short-circuit: do NOT continue through QC with a bad background
            generate_manifest(config, ValidationStatus.FAIL,
                              background_file_result.issues, None, None,
                              args.output, start_time, time.time(),
                              input_path=args.input,
                              normalized_path=norm_result.normalized_path,
                              background_path=args.background)
            print_failure(background_file_result)
            sys.exit(1)
        background_norm = bg_adapter.normalize(args.background, args.output,
                                                label="background")
    
    # ── Step 8: Validate + prepare custom gene sets ───────────────────
    # Called UNCONDITIONALLY so the legal-combinations table is enforced:
    #   --gene-sets custom + no files  → FAIL
    #   --custom-gene-set-action without files → FAIL
    #   default + no custom → PASS (no-op)
    # Returns GeneSetPreparationResult with normalized_paths — the actual
    # files the engine should receive (converted from GMT if needed).
    custom_gs_prep = prepare_custom_gene_sets(
        custom_files=args.custom_gene_set_files or [],
        custom_format=config.get("custom_gene_set_format"),
        custom_action=config.get("custom_gene_set_action"),
        profile=config["gene_sets"],
        output_dir=args.output)
    # Short-circuit: do NOT proceed to path resolution with bad config
    if custom_gs_prep.status == ValidationStatus.FAIL:
        generate_manifest(config, ValidationStatus.FAIL,
                          custom_gs_prep.issues, None, None, args.output,
                          start_time, time.time(),
                          args.input, norm_result.normalized_path)
        print_failure(custom_gs_prep)
        sys.exit(1)
    
    # ── Step 9: Resolve final gene-set + reference paths ────────────
    # Uses custom_gs_prep.normalized_paths (the converted files),
    # NOT args.custom_gene_set_files (the originals).
    # custom_gs_prep is always set (called unconditionally above).
    base_dir = resolve_base_dir()
    ref_paths = resolve_reference_paths(config["genome_build"], base_dir)
    gs_paths = resolve_gene_set_paths(
        config["gene_sets"],
        custom_gs_prep.normalized_paths,
        config.get("custom_gene_set_action"),
        base_dir)
    gene_map_path = resolve_gene_map_path(base_dir)
    
    # ── Step 10: Gene / reference QC on NORMALIZED genes ────────────
    gene_map = load_gene_map(gene_map_path)
    gene_locs = load_gene_locations(ref_paths["gene_loc"])
    gene_qc = validate_normalized_genes(norm_result, gene_map, gene_locs, gs_paths)
    write_gene_qc_tsv(gene_qc, args.output)
    
    # ── Step 11: Background QC (consistency with input genes) ───────
    if background_norm:
        background_qc = validate_background(background_norm, norm_result, gene_map)
    
    # ── Step 12: Combined validation status ─────────────────────────
    combined_status, all_issues = compute_combined_status(
        file_result, norm_result, gene_qc,
        background_file_result=background_file_result,
        background_qc_result=background_qc,
        custom_gs_prep_result=custom_gs_prep)
    if combined_status == ValidationStatus.FAIL:
        generate_manifest(config, combined_status, all_issues, gene_qc, None,
                          args.output, start_time, time.time(),
                          input_path=args.input,
                          normalized_path=norm_result.normalized_path,
                          gene_map_path=gene_map_path,
                          gene_set_paths=gs_paths,
                          reference_paths=ref_paths,
                          background_qc=background_qc,
                          background_path=args.background,
                          background_norm_path=(
                              background_norm.normalized_path
                              if background_norm else None))
        print_failure_with_qc(combined_status, all_issues, args.output)
        sys.exit(1)
    
    # ── Step 13–15: Build + run engine ─────────────────────────────
    cmd = build_priors_command(config, adapter, norm_result, ref_paths, gs_paths,
                               gene_map_path, args.output, background_norm)
    exit_code = run_engine(cmd, args.output)
    
    # ── Step 16–17: Final manifest + verified report ────────────────
    end_time = time.time()
    generate_manifest(config, combined_status, all_issues, gene_qc, exit_code,
                      args.output, start_time, end_time,
                      input_path=args.input,
                      normalized_path=norm_result.normalized_path,
                      gene_map_path=gene_map_path,
                      gene_set_paths=gs_paths,
                      reference_paths=ref_paths,
                      background_qc=background_qc,
                      background_path=args.background,
                      background_norm_path=(
                          background_norm.normalized_path
                          if background_norm else None))
    if exit_code == 0:
        generate_text_report(config, gene_qc, ..., args.output)
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
```

#### Step 9: Dockerfile Update

```dockerfile
FROM ubuntu:plucky
RUN apt -y update && apt -y upgrade && apt -y install python3 python3-scipy
WORKDIR /app
# Reference data
ADD data/gene_set_list_mouse_2024.txt data/
ADD data/gene_set_list_msigdb_nohp.txt data/
ADD data/portal_gencode.gene.map data/
ADD data/NCBI37.3.plink.gene.loc data/
ADD data/refGene_hg19_TSS.subset.loc data/
ADD data/NCBI37.3.plink.gene.exons.loc data/
# Engine (unchanged)
ADD priors.py priors.py
# Wrapper (new)
ADD run_pigean.py run_pigean.py
ADD pigean/ pigean/

ENTRYPOINT ["python3", "-u", "/app/run_pigean.py"]
```

> **Note:** No extra packages needed — using JSON (stdlib) for all config/manifest files.

#### Step 10: WDL Update

```wdl
version 1.0

workflow rock_pigean {
    input {
        File input_file
        String analysis_type = "positive-controls"
        String preset = "standard"
        String gene_set_profile = "default"
        String genome_build = "hg19"
        Array[File] custom_gene_sets = []
        String custom_gene_set_format = "engine"
        String custom_gene_set_action = "supplement"
        File? background_gene_list
        Int memory_gb = 8
        Int cpu = 4
        Int disk_gb = 50
    }
    call run_pigean { input:
        input_file = input_file,
        analysis_type = analysis_type,
        preset = preset,
        gene_set_profile = gene_set_profile,
        genome_build = genome_build,
        custom_gene_sets = custom_gene_sets,
        custom_gene_set_format = custom_gene_set_format,
        custom_gene_set_action = custom_gene_set_action,
        background_gene_list = background_gene_list,
        memory_gb = memory_gb,
        cpu = cpu,
        disk_gb = disk_gb,
    }
    output {
        File gs = run_pigean.gs
        File gss = run_pigean.gss
        File ggss = run_pigean.ggss
        File params = run_pigean.params
        File manifest = run_pigean.manifest
        File gene_qc = run_pigean.gene_qc
        File report = run_pigean.report
        File resolved_config = run_pigean.resolved_config
    }
}

task run_pigean {
    input {
        File input_file
        String analysis_type
        String preset
        String gene_set_profile
        String genome_build
        Array[File] custom_gene_sets
        String custom_gene_set_format
        String custom_gene_set_action
        File? background_gene_list
        Int memory_gb
        Int cpu
        Int disk_gb
    }
    runtime {
        docker: "gcr.io/nitrogenase-docker/rock-pigean:2.0.0"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
    command <<<
        python3 -u /app/run_pigean.py \
            --analysis ~{analysis_type} \
            --input ~{input_file} \
            --gene-sets ~{gene_set_profile} \
            --genome-build ~{genome_build} \
            --preset ~{preset} \
            --output results \
            ~{if length(custom_gene_sets) > 0 then "--custom-gene-set-files " + sep(" ", custom_gene_sets) + " --custom-gene-set-format " + custom_gene_set_format + " --custom-gene-set-action " + custom_gene_set_action else ""} \
            ~{if defined(background_gene_list) then "--background " + background_gene_list else ""}
    >>>
    output {
        File gs = "results/gs.out"
        File gss = "results/gss.out"
        File ggss = "results/ggss.out"
        File params = "results/p.out"
        File manifest = "results/run_manifest.json"
        File gene_qc = "results/input_gene_qc.tsv"
        File report = "results/report.txt"
        File resolved_config = "results/resolved_config.json"
    }
}
```

---

## Phase 2 — Custom Gene Sets + Background Genes (~1.5 days)

### Output Structure

```
outputs/phase2_gene_sets/
├── default/                          ← baseline (should match golden)
│   ├── input/
│   ├── normalized/
│   ├── priors_command.txt
│   ├── resolved_config.json
│   ├── gs.out, gss.out, ggss.out, p.out
│   └── report.txt
│
├── msigdb_only/
│   └── ...
│
├── mouse_only/
│   └── ...
│
├── custom_engine_format/
│   ├── input/
│   │   └── custom_pathways.txt
│   ├── normalized/
│   │   └── custom_pathways.txt       ← same (no conversion needed)
│   └── ...
│
├── custom_gmt/
│   ├── input/
│   │   └── custom_pathways.gmt       ← original GMT format
│   ├── normalized/
│   │   └── custom_pathways.engine.txt ← converted (description column stripped)
│   └── ...
│
├── custom_supplement/                ← custom pathways ADDED to defaults
│   └── ...
│
├── custom_replace/                   ← custom pathways REPLACING defaults
│   └── ...
│
├── background_test/
│   ├── without_background/
│   ├── with_background/
│   │   ├── input/
│   │   │   ├── positive_controls.txt
│   │   │   └── background_original.txt
│   │   ├── normalized/
│   │   │   ├── positive_controls.txt
│   │   │   └── background.txt
│   │   └── ...
│   └── comparison/
│       ├── gene_score_comparison.tsv
│       └── enrichment_change_summary.txt
│
└── summary.txt
```

### Step 1: `pigean/adapters/gene_sets.py`

> **GMT format is declared explicitly, not auto-detected.** v1 requires `--custom-gene-set-format`.

> **Replace vs supplement is declared explicitly.** `--custom-gene-set-action replace` uses only the custom files; `--custom-gene-set-action supplement` adds them alongside the **base profile** (whichever `--gene-sets` selected). Recorded in `resolved_config.json`, `run_manifest.json`, and `report.txt`.

#### Legal Gene-Set Combinations

There is **one unambiguous rule**: `--gene-sets` selects the base profile, `--custom-gene-set-action` decides whether custom files replace or supplement that base.

| `--gene-sets` | `--custom-gene-set-files` | `--custom-gene-set-action` | Result |
|---|---|---|---|
| `default` | *(none)* | *(none)* | Mouse + MSigDB |
| `mouse-only` | *(none)* | *(none)* | Mouse only |
| `msigdb-only` | *(none)* | *(none)* | MSigDB only |
| `custom` | required | `replace` | Custom only |
| `custom` | required | `supplement` | Default (Mouse + MSigDB) + Custom |
| `mouse-only` | provided | `supplement` | Mouse + Custom |
| `msigdb-only` | provided | `supplement` | MSigDB + Custom |
| `default` | provided | `replace` | Custom only (base ignored) |
| *(any)* | provided | *(missing)* | **ERROR** — action is required |
| *(any)* | *(none)* | provided | **ERROR** — action without files is meaningless |
| `custom` | *(none)* | *(any)* | **ERROR** — custom profile requires files |

These rules are **unit-tested in `test_gene_sets.py`**.

> **WDL default:** The WDL sets `custom_gene_set_action = "supplement"` so that submitting a workflow with custom files but no explicit action choice adds them alongside defaults (the least surprising behavior for biologists). This is a deliberate convenience default in the WDL layer, not a silent assumption — the CLI still requires the flag explicitly.

```python
# pigean/adapters/gene_sets.py
#
# Responsibility: validate and convert custom gene-set files.
# Does NOT assemble the final --X-in list — that is
# references.py::resolve_gene_set_paths().

def convert_gmt_to_engine_format(gmt_path, output_path):
    """Convert GMT (name→desc→genes) to engine format (name→genes) by stripping
    column index 1 (the description field). Validates that no description strings
    ended up in the gene list by checking against the gene map after conversion."""
    ...
```

> **Note:** `prepare_custom_gene_sets()` lives in `validation.py` (it orchestrates validation + conversion). `convert_gmt_to_engine_format()` lives in `gene_sets.py` (the low-level conversion). `resolve_gene_set_paths()` lives in `references.py` (assembles base + custom into the final `--X-in` list). One function per job, no overlap.

### Step 2: Background Gene Support

Background files go through the **same validate → normalize pipeline** as input genes:

```
args.background
      ↓
validate_file()
      ↓
normalize() → input/background_original.txt + normalized/background.txt
      ↓
validate_background() → check consistency with input genes
      ↓
engine gets: --positive-controls-all-in <normalized/background.txt>
```

The engine receives the **normalized** background path, not the raw config value. The `input/` and `normalized/` directories show both the original and cleaned versions for audit.

### Step 3: Update Validation

- Validate custom gene set files exist and are readable
- Require `--custom-gene-set-format` when `--custom-gene-set-files` is used (fail otherwise)
- Require `--custom-gene-set-action` when `--custom-gene-set-files` is used (fail otherwise)
- Warn if custom gene sets have very few entries
- Validate background file genes against gene map
- **Report input genes missing from declared background** (warn, do not fail — policy TBD after confirming engine behavior)
- Record custom gene set action and background coverage in manifest and report

---

## Dependencies

**No new dependencies for v1.** Using JSON (stdlib) instead of YAML for all configuration and manifest files. Everything the wrapper needs:

| Package | Source | Used For |
|---------|--------|----------|
| `json` | stdlib | Manifest, resolved config, checkpoints |
| `argparse` | stdlib | CLI parsing |
| `subprocess` | stdlib | Engine execution |
| `hashlib` | stdlib | File checksums |
| `dataclasses` | stdlib (3.7+) | Config/result types |
| `shutil` | stdlib | Copy input files |
| `platform` | stdlib | Python version in manifest |
| `time` | stdlib | Run timing |
| `scipy` | engine dep | (Spearman correlation in golden test comparison scripts) |
| `numpy` | engine dep | Engine dependency |

---

## Verification Plan

### After Phase 0
- `run_golden_test.sh` passes using **tolerance-based comparison**:
  - Repeated runs determine which quantities are deterministic vs stochastic
  - Deterministic quantities → exact checks
  - Stochastic quantities → empirical tolerances from observed variability
  - Ranking correlation → Spearman ρ threshold
- 5 repeated runs in `outputs/phase0_golden/` with `baseline_summary/` computed
- Tolerances in `checkpoints.json` calibrated from empirical variability
- `num_gibbs_iter` in `p.out` clarified: is 499 a zero-based index or the count?

### After Phase 1
- **Parameter comparison:** `outputs/phase1_wrapper/comparison/parameter_comparison.txt` confirms the wrapper produces the same effective engine parameters and reference files as the original command (allowing expected path differences in `command_diff.txt`)
- **Result compatibility:** `outputs/phase1_wrapper/comparison/compatibility_result.txt` = PASS (using the same tolerance-based checks as the golden test)
- **Manifest test:** `run_manifest.json` is valid JSON with all required fields including checksums, timing, and combined validation status
- **Failed run audit trail:** Running with a nonexistent input path → `resolved_config.json` + manifest still written with `execution_status: "NOT_RUN"`, clear error
- **Requested-but-invalid background:** Running with `--background /nonexistent` → combined status = `FAIL`, not silently ignored
- **Validation test:** Running with a file containing only `FAKEGENE1` → combined status = `FAIL`, manifest says `NOT_RUN`, QC file written
- **QC test:** `input_gene_qc.tsv` correctly identifies PLINK1 as unresolved, describes **normalized** genes
- **Normalized input test:** `normalized/` contains the cleaned input, `input/` contains the original
- **Report test:** Text report displays only metrics with verified interpretation contract entries
- **Non-empty output dir test:** Running into an existing output directory → fails unless `--overwrite`; with `--overwrite`, old contents are deleted before the new run starts (no stale files survive)
- **Config precedence test:** `test_config.py` proves DEFAULTS < config JSON < explicit CLI
- **Engine command test:** `test_engine_command.py` proves wrapper generates correct effective parameters

### After Phase 2
- `outputs/phase2_gene_sets/default/` passes the golden test (tolerance-based)
- `outputs/phase2_gene_sets/custom_gmt/normalized/` has description column stripped
- Using `--custom-gene-set-files` without `--custom-gene-set-format` → clear error
- Using `--custom-gene-set-files` without `--custom-gene-set-action` → clear error
- `custom_supplement/` includes both default and custom gene sets; `custom_replace/` includes only custom
- **Gene-set combination rules:** `test_gene_sets.py` covers the full legal-combinations table — including `msigdb-only + supplement + custom → MSigDB + Custom`, `--custom-gene-set-files` without `--custom-gene-set-action` → error, `--gene-sets custom` without files → error
- `resolved_config.json` records the `custom_gene_set_action` value
- `outputs/phase2_gene_sets/background_test/comparison/` successfully **quantifies changes** in gene/pathway statistics with versus without the declared background (any direction of change is a valid finding)
- Background QC reports any input genes missing from declared background
- Background file appears in `input/` (original) and `normalized/` (cleaned)
- WDL correctly passes custom gene sets, format, action, and background file to `run_pigean.py`
- Golden test still passes with default settings

### End-to-end
```bash
# Should produce outputs with full audit trail
python3 run_pigean.py \
    --analysis positive-controls \
    --input ex/gene_list \
    --output /tmp/v1_test

# Verify output structure
ls /tmp/v1_test/
# → input/  normalized/  priors_command.txt  resolved_config.json
#   input_gene_qc.tsv  run_manifest.json  pigean_run.log
#   gs.out  gss.out  ggss.out  p.out  report.txt

# Validate against golden using tolerance-based checks (not raw diff)
bash tests/golden/run_golden_test.sh /tmp/v1_test
# → PASS (structural counts exact, gene priors within tolerance,
#         top-N overlap above threshold, rank correlation above threshold)
```

---

## .gitignore Additions

```
.idea/
outputs/
runs/
```

The `outputs/` and `runs/` directories contain generated artifacts — never checked into version control. The frozen golden reference files in `tests/golden/expected/` *are* checked in.
