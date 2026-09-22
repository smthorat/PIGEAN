# PIGEAN Pipeline — Generalization Plan

> **Design principle:** A generalized PIGEAN is not "expose every `priors.py` flag in WDL."
> It is: **a user can bring a supported type of biological evidence, select reference data, run one simple interface, and receive validated, interpretable results.**

> **Revision note (2026-09-10, rev 2):** Incorporates two rounds of design review feedback. Key changes: (1) basic report and reproducibility moved into v1; (2) `optparse` corrected from `argparse`; (3) "discovered genes" renamed to "additional prioritized candidates"; (4) background gene framing corrected — genes outside the positive list are not known negatives; (5) `--correct-huge` reclassified as needing verification; (6) adapter contract expanded to return mapping decisions and exclusions, with explicit ambiguous-alias policy; (7) release order revised to v1/v1.1/v2/v2.1/later; (8) input gene default probability (≈0.95) documented; (9) GMT format support decided; (10) stochastic tolerance estimation from repeated runs required before fixing thresholds; (11) probability separated from BF scales — different conversion with background-odds dependency; (12) engine-freeze rule qualified — wrapper changes to filtering/mapping/references require separate validation; (13) manifest expanded with dependency versions, checksums, three-status model, exit code, and log locations; (14) timeline corrected to 37–43 working days (7.5–9 weeks); (15) shell examples fixed for copy-paste safety.

---

## Table of Contents

1. [Current State](#1-current-state)
2. [Target Architecture](#2-target-architecture)
3. [What priors.py Already Supports (Hidden Capabilities)](#3-what-priorspy-already-supports)
4. [Phased Execution Plan](#4-phased-execution-plan)
5. [Phase 0 — Golden Test and Specification](#5-phase-0--golden-test-and-specification)
6. [Phase 1 — Wrapper Layer + Parameterized WDL](#6-phase-1--wrapper-layer--parameterized-wdl)
7. [Phase 2 — Custom Gene Sets + Background Genes](#7-phase-2--custom-gene-sets--background-genes)
8. [Phase 3 — Multi-Evidence Input Support](#8-phase-3--multi-evidence-input-support)
9. [Phase 4 — GWAS + Genome Build Reference Packs](#9-phase-4--gwas--genome-build-reference-packs)
10. [Phase 5 — Full Interpretation Report and Convergence Diagnostics](#10-phase-5--full-interpretation-report-and-convergence-diagnostics)
11. [Phase 6 — Factor / PheWAS / Advanced Modes](#11-phase-6--factor--phewas--advanced-modes)
12. [Phase 7 — Production Infrastructure](#12-phase-7--production-infrastructure)
13. [Reference Pack Strategy](#13-reference-pack-strategy)
14. [Evidence Adapter Architecture](#14-evidence-adapter-architecture)
15. [User Interface Design — Normal vs. Advanced](#15-user-interface-design--normal-vs-advanced)
16. [Known Warnings — Verified vs. Hypothesized](#16-known-warnings--verified-vs-hypothesized)
17. [Effort Estimates — Implementation vs. Validation](#17-effort-estimates--implementation-vs-validation)
18. [MVP Definitions](#18-mvp-definitions)
19. [Open Questions](#19-open-questions)
20. [How to Run — Manual Commands for Each Phase](#how-to-run--manual-commands-for-each-phase)

---

## 1. Current State

### What We Have Today

```
User provides:   gene_list (File), output_files_base_name (String)
                        │
                        ▼
              ┌─────────────────┐
              │  rock_pigean.wdl │  ← Only 2 inputs exposed
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │  Docker image    │  ← All reference data baked in
              │  priors.py gibbs │  ← Only 1 mode used (of 11)
              │  + 6 data files  │  ← Only hg19, only 2 gene set DBs
              └────────┬────────┘
                       │
                       ▼
              4 output files: gs, gss, ggss, p
              (no QC, no report, no manifest)
```

### What Is Hardcoded

| Item | Hardcoded Value | Impact |
|---|---|---|
| PIGEAN mode | `gibbs` | Only 1 of 11 modes accessible |
| Input type | `--positive-controls-in` only | No GWAS, exome, or gene score input |
| Gene set databases | 2 files baked into Docker | No custom gene sets without rebuild |
| Genome build | hg19/GRCh37 only | No hg38 support |
| Gene map | `portal_gencode.gene.map` | No alternative annotations |
| Max gene sets | `5000` | Not tunable |
| Filter thresholds | `gene=1, gene_set=0.01` | Not tunable |
| Debug level | `3` | Not tunable |
| Runtime resources | Not declared | Jobs may OOM silently |
| Docker base image | `ubuntu:plucky` (mutable release tag, not pinned digest) | Not strictly reproducible |
| Background genes | Not provided | Enrichment may be inflated (unverified — see §16) |
| Input validation | None | Bad genes silently score zero |
| Run provenance | None | No record of software/reference versions |

---

## 2. Target Architecture

### End-State Design

```
                         USER
                           │
                 ┌─────────▼──────────┐
                 │  PIGEAN interface   │  ← Simple parameters
                 │  (run_pigean.py)    │  ← Analysis types, not engine flags
                 └─────────┬──────────┘
                           │
           ┌───────────────┼────────────────┐
           │               │                │
     positive controls    GWAS            exome / gene scores
           │               │                │
           └────── evidence adapters ───────┘
                           │
                           ▼
                    Gene-level evidence
                    (normalized / engine-ready)
                           │
          ┌────────────────┼────────────────┐
          │                                 │
   Default gene sets                 Custom gene sets
   (mouse + MSigDB)                  (user-provided)
          │                                 │
          └────────────────┬────────────────┘
                           │
                 ┌─────────▼──────────┐
                 │   Reference pack   │  ← Versioned, with manifest
                 │   (hg19 or hg38)   │
                 └─────────┬──────────┘
                           │
                 ┌─────────▼──────────┐
                 │   PIGEAN engine    │  ← priors.py (untouched)
                 │   (Gibbs sampler)  │
                 └─────────┬──────────┘
                           │
          ┌────────────────┼─────────────────┐
          ▼                ▼                  ▼
    gs/gss/ggss/p    QC + manifest      HTML report
    (raw outputs)    (run_manifest.json) (interpretation)
```

### Why a Wrapper

Currently:
```
WDL → priors.py + dozens of flags
```

Proposed:
```
WDL / Nextflow / Docker / CLI / Web
                │
                ▼
          run_pigean.py        ← User-facing layer
                │
                ▼
          priors.py gibbs      ← Engine (unchanged)
          + internal flags
```

The wrapper gives us:

1. **One interface, many execution platforms.** WDL, Nextflow, Docker, Slurm, or a future web app all call the same `run_pigean.py`.
2. **User-facing analysis types** instead of engine flags. Scientists choose `--analysis positive-controls`, not `--positive-controls-in` plus 15 flags.
3. **Input validation before expensive compute.** Catch bad gene names, missing columns, mismatched builds before a 30-minute Gibbs run.
4. **Run manifests and reports automatically.** Every run records what was done and produces interpretable output.
5. **Validated defaults.** Users don't need to decide between 7 or 10 chains.

---

## 3. What priors.py Already Supports

The engine (`priors.py`, 16,568 lines) already has extensive capabilities that the wrapper will expose incrementally.

### Modes (11 available, 1 currently used)

| Mode | Internal Use | User-Facing Preset |
|---|---|---|
| `gibbs` | Full Bayesian Gibbs sampling | `--preset standard` |
| `naive_priors` | Fast approximate priors | `--preset fast` (**experimental** — see §16) |
| `factor` / `naive_factor` | Pathway decomposition | `--analysis factor` |
| `priors` | Full pipeline, no sampling | Advanced only |
| `beta_tildes` | Marginal associations only | Advanced only |
| `betas` | Conditional betas only | Advanced only |
| `sigma` | Learn hyperparameter | Advanced only |
| `sim` | Simulations | Advanced only |
| `pops` / `naive_pops` | PoPS-style | Advanced only |

### Input Evidence Types (multiple supported, 1 currently used)

The table below lists known input flags. Exact counts should be confirmed by programmatic inventory of `priors.py`'s `optparse` definitions (the engine uses `optparse`, not `argparse`).

| Evidence Type | priors.py Flag | Adapter Complexity |
|---|---|---|
| Gene list | `--positive-controls-in` | Simple (in use) |
| Gene list (inline) | `--positive-controls-list` | Simple |
| Background genes | `--positive-controls-all-in` | Simple |
| Gene Bayes factors | `--gene-bfs-in` | Medium — must distinguish raw BF, natural-log BF, log10 BF, and probability inputs; default reader expects `Gene` and `log_bf` columns |
| Gene Z-scores | `--gene-zs-in` | Medium — reader contains multiple interpretations (raw log-odds, probability mapping); validate selected conversion before supporting arbitrary Z-scores |
| Gene percentiles | `--gene-percentiles-in` | Simple |
| Exome associations | `--exomes-in` | Medium |
| GWAS summary stats | `--gwas-in` | **Complex** (coordinates, builds, S2G, LD) |
| Credible sets | `--credible-sets-in` | Complex |
| SNP-to-gene map | `--s2g-in` | Complex |

### Output Types (multiple available, 4 currently used)

The 4 core outputs (gs, gss, ggss, p) plus additional diagnostic and advanced outputs: trace files (convergence diagnostics), overlap stats, gene effectors, pheWAS stats, factors, clusters, gene-pheno stats, X matrix, among others. Exact counts should be confirmed by `optparse` inventory.

---

## 4. Phased Execution Plan

```
Phase 0 ─── Golden test + specification
    │
    ▼
Phase 1 ─── Wrapper (run_pigean.py) + parameterized WDL + input validation
    │
    ▼
Phase 2 ─── Custom gene sets + background genes
    │
    ▼
Phase 3 ─── Multi-evidence input support (BF / Z / exome)
    │
    ▼
Phase 4 ─── GWAS + hg19/hg38 reference packs
    │
    ▼
Phase 5 ─── Automatic interpretation / HTML report
    │
    ▼
Phase 6 ─── Factor / PheWAS / advanced modes
    │
    ▼
Phase 7 ─── Production infrastructure
```

### Why This Order

- **Phase 0 first** because every subsequent change needs a regression test.
- **Wrapper (Phase 1) before features** because it establishes the architecture everything else plugs into.
- **Custom gene sets (Phase 2) before GWAS (Phase 4)** because gene sets are the annotation layer that all analyses share, and it's simpler.
- **Multi-evidence inputs (Phase 3) before GWAS (Phase 4)** because GWAS is the most complex evidence adapter (coordinates, builds, LD, S2G mapping). Gene BFs and Z-scores are much simpler and prove the evidence-adapter pattern works. Then GWAS becomes "just another adapter."
- **Report (Phase 5) before advanced modes (Phase 6)** because interpretability benefits all users, while advanced modes benefit few.

---

## 5. Phase 0 — Golden Test and Specification

### Goal

Freeze the current known-good run as a regression test. Every future change must prove:

> **Generalization changed usability, not scientific behavior.**

### Deliverables

```
tests/
├── golden/
│   ├── input/
│   │   └── positive_controls.txt        # The 10-gene list
│   ├── expected/
│   │   ├── gs.golden.out
│   │   ├── gss.golden.out
│   │   ├── ggss.golden.out
│   │   └── p.golden.out
│   ├── checkpoints.yaml                 # Key values to validate
│   └── run_golden_test.sh               # Runs pipeline, compares
```

> **Golden test does NOT require byte-for-byte equality for stochastic outputs.**
> Gibbs sampling is inherently random. The engine calls `random.seed(0)`, but NumPy random sampling is **not** correspondingly seeded in the inspected file — this does **not** establish reproducible Gibbs sampling. Unless a full reproducible seed pathway is confirmed during Phase 0, the golden test validates structural invariants and statistical tolerances — not exact file contents.

If a complete seed pathway is confirmed, record and fix it for the golden run. If not, validate:
- Row counts (exact — verify which counts are truly deterministic via repeated runs)
- Specific gene priors within tolerance
- Top-N gene overlap (e.g., top 5 candidate genes by prior)
- Top-N pathway overlap
- Rank correlation of gene priors (target Spearman ρ ≥ 0.99, but **estimate actual stochastic tolerances from repeated baseline runs before fixing thresholds**)
- Model hyperparameters within tolerance

> **Important:** Tolerance thresholds (Spearman ≥ 0.99, ±0.05, etc.) listed here and in `checkpoints.yaml` are provisional. Run the golden test case multiple times to measure empirical variability before fixing acceptance criteria.

### Three Types of Validation Checks

Use the golden run to establish three layers of confidence:

| Check | What It Establishes |
|---|---|
| **Compatibility:** Original `priors.py` command vs. wrapper, with identical inputs and references | The wrapper preserves existing behavior |
| **Variation:** Small, medium, partially-mapped, and custom-annotation examples | The workflow handles supported input variation |
| **Stability:** Repeated runs, appropriate null controls, and held-out evidence | Rankings are stable and scientific performance is credible |

> **Circularity warning:** For prediction benchmarks using held-out evidence, verify that the held-out evidence has not already entered the gene-set annotations. Otherwise, apparently successful recovery can be circular.

### Checkpoint Values to Track

```yaml
# checkpoints.yaml — values that must remain stable within tolerance
# NOTE: All tolerances below are provisional estimates.
# Run the golden case ≥5 times to measure empirical variability before fixing these.

gene_priors:
  LEP:     {expected: 1.89,   tolerance: 0.05}
  ADIPOQ:  {expected: 1.36,   tolerance: 0.05}
  BRD2:    {expected: 0.803,  tolerance: 0.05}
  LECT2:   {expected: 0.197,  tolerance: 0.05}
  PLINK1:  {expected: -0.010, tolerance: 0.02}

candidate_genes:
  LEPR:    {expected: 0.988,  tolerance: 0.05}
  PPARG:   {expected: 0.640,  tolerance: 0.05}

model:
  num_genes:               {expected: 42216, tolerance: 0}
  num_gene_sets_retained:  {expected: 2094,  tolerance: 10}
  p_version1:              {expected: 0.00356, tolerance: 0.001}
  sigma2_version1:         {expected: 6.41e-08, tolerance: 1e-08}

counts:
  gs_rows:    {expected: 42216, tolerance: 0}
  gss_rows:   {expected: 41627, tolerance: 0}   # structural — verify determinism via repeated runs
```

### Why This Matters

We already discovered during our initial runs that some output interpretations were not immediately obvious (e.g., what `combined_D` means, how `prior` vs `combined` differ). When modifying the wrapper, the golden test proves that scientific behavior is preserved.

### Engine-Freeze Rule

> **`priors.py` is not modified during v1 (Phases 0–2).**

All generalization work happens in the wrapper layer (`run_pigean.py`), WDL, Docker, and supporting scripts. The engine is treated as a fixed dependency. This ensures:

1. The golden test baseline remains stable — the engine itself is unchanged. **Compatibility is assessed with identical evidence, references, preprocessing, and settings. Changes to those inputs require separate validation.**
2. Bugs found in the wrapper are clearly wrapper bugs, not engine regressions.
3. If we later need to patch `priors.py`, the golden test catches any scientific drift.

Engine changes (if needed) are gated on a separate PR with its own before/after golden-test comparison.

**Estimated effort:** 2–3 hours implementation + 0.5 day for repeated baseline runs to establish stochastic tolerances.

---

## 6. Phase 1 — Wrapper Layer + Parameterized WDL

### Goal

Build `run_pigean.py` — a user-facing interface that translates simple parameters into `priors.py` engine flags. Move input validation here (not Phase 6). Parameterize the WDL to call the wrapper.

### run_pigean.py Interface

```bash
# Simple usage — matches current behavior exactly
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --gene-sets default \
    --genome-build hg19 \
    --output results/run1

# Internally translates to:
priors.py gibbs \
    --X-in /app/data/gene_set_list_mouse_2024.txt \
    --X-in /app/data/gene_set_list_msigdb_nohp.txt \
    --gene-map-in /app/data/portal_gencode.gene.map \
    --positive-controls-in genes.txt \
    --max-num-gene-sets 5000 \
    --gene-filter-value 1 \
    --gene-set-filter-value 0.01 \
    --gene-loc-file /app/references/hg19/gene.loc \
    --gene-loc-file-huge /app/references/hg19/tss.loc \
    --exons-loc-file-huge /app/references/hg19/exons.loc \
    --gene-stats-out results/run1/gs.out \
    --gene-set-stats-out results/run1/gss.out \
    --gene-gene-set-stats-out results/run1/ggss.out \
    --params-out results/run1/p.out \
    --debug-level 3
```

### User-Facing Analysis Types (Not Raw Modes)

| User selects | Internally runs | Description |
|---|---|---|
| `--analysis positive-controls` | `gibbs --positive-controls-in` | Enrich from a gene list |
| `--analysis gwas` | `gibbs --gwas-in` | Enrich from GWAS summary stats (Phase 4) |
| `--analysis exome` | `gibbs --exomes-in` | Enrich from exome associations (Phase 3) |
| `--analysis gene-scores` | `gibbs --gene-bfs-in` / `--gene-zs-in` | Enrich from precomputed scores (Phase 3) |
| `--analysis factor` | `factor` mode | Pathway decomposition (Phase 6) |

### User-Facing Presets (Not Raw Tuning Parameters)

| User selects | Internally sets | Description | Status |
|---|---|---|---|
| `--preset standard` (default) | `gibbs`, current defaults | Current validated pipeline configuration | **Supported** |
| `--preset fast` | `naive_priors` | Approximate alternative | **Experimental** — pending §16 validation |
| `--preset thorough` | `gibbs`, TBD settings | Extended sampling configuration | **Experimental** — settings TBD pending benchmarks |

> **v1 ships with `standard` only.** The `fast` and `thorough` presets are not exposed in v1. They become available (gated behind `--experimental`) once the benchmarks in §16 are complete.

Scientists should not need to decide between 7 or 10 chains. Validated defaults handle that.

### Advanced Parameters (Separate Namespace)

For power users only, via a `--config` YAML file:

```bash
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --config advanced.yaml
```

The wrapper merges the config with validated defaults, writes `resolved_config.yaml` to the output directory, and translates to engine flags. These parameters are **not** part of the normal interface.

### Input Validation (Day 1, Not Phase 6)

The wrapper validates inputs **before** calling `priors.py`:

```
INPUT VALIDATION
────────────────
Input genes:                 10
Recognized gene symbols:      <value>  (found in gene map)
With genomic coordinates:     <value>  (mapped to hg19 loci)
With ≥1 gene-set membership:  <value>  (have pathway annotations)
Unresolved genes:             1 ⚠
  PLINK1 — unresolved in one or more reference layers;
           exact failure reported by input_gene_qc.tsv

Gene sets retained for modeling: 2,094
Genome build:             hg19
Validation outcome:       PASS_WITH_WARNINGS ⚠

Proceeding with <value>/10 resolved genes.
```

For GWAS (Phase 4):
```
GWAS INPUT VALIDATION
─────────────────────
Genome build:             hg38
Variants:             1,247,893
Columns detected:     CHR, BP, P, BETA, SE, N
Missing beta:                 0 ✓
Missing SE:                  12 ⚠ (dropped)
Chromosome format:       1-22,X ✓
Coordinate range:     consistent with hg38 ✓
```

This is essential scientific QC, not an optional production feature.

#### Validation Outcomes

The wrapper assigns a formal validation outcome, written to the run manifest and printed at the start of every run:

| Outcome | Meaning | Action |
|---|---|---|
| `PASS` | All input genes resolved, all QC checks within normal range | Run proceeds |
| `PASS_WITH_WARNINGS` | Run can proceed but some inputs degraded (e.g., unmapped genes, small gene list) | Run proceeds; warnings in manifest and report |
| `FAIL` | Run cannot produce meaningful results (e.g., 0 mapped genes, no gene sets loaded, unrecognized analysis type) | Run aborts with clear error message |

v1 hard-fail criteria (clearly unusable):
- 0 recognized genes → `FAIL`
- No usable gene sets loaded → `FAIL`
- Invalid input format (e.g., binary file, missing required columns) → `FAIL`
- Incompatible requested references (e.g., unknown genome build) → `FAIL`

For partial mapping (some genes unresolved but at least 1 usable), issue `PASS_WITH_WARNINGS` and record the mapped fraction. A percentage-based failure threshold (e.g., <50%) can be introduced after benchmarking different gene-list sizes establishes a sensible cutoff.

#### Input Gene QC Output

Every run writes `input_gene_qc.tsv` to the output directory — a per-gene record of how each input gene was resolved across the reference layers:

```
# Illustrative; recognized/coordinate status will be determined by Phase 1 validator.
input_gene	recognized	coordinates	annotation_memberships	status
LEP	yes	yes	930	PASS
ADIPOQ	yes	yes	650	PASS
BRD2	yes	yes	432	PASS
SOX2	yes	yes	499	PASS
LITAF	yes	yes	453	PASS
LECT2	yes	yes	114	PASS
SLCO1B1	yes	yes	119	PASS
RMI2	yes	yes	164	PASS
CD300LG	yes	yes	66	PASS
PLINK1	<value>	<value>	0	UNRESOLVED
```

This file is always written. It distinguishes which lookup step failed for each gene (gene-map recognition, coordinate mapping, gene-set membership), rather than collapsing everything into a binary "mapped/unmapped." Downstream tools or reports can consume it without parsing log output.

### Run Manifest (Day 1)

Every run produces `run_manifest.json`:

```json
{
  "wrapper_version": "1.0.0",
  "engine": {
    "name": "priors.py",
    "version": "<git commit SHA or file SHA-256>"
  },
  "dependencies": {
    "python": "3.11.9",
    "numpy": "1.26.4",
    "scipy": "1.13.0"
  },
  "container_digest": "sha256:abc123...",
  "analysis_type": "positive_controls",
  "preset": "standard",
  "genome_build": "hg19",
  "reference_pack": "PIGEAN-GRCh37-v1",
  "reference_checksums": {
    "gene.loc": "sha256:...",
    "tss.loc": "sha256:...",
    "exons.loc": "sha256:..."
  },
  "gene_sets": ["mouse_2024", "msigdb_nohp"],
  "gene_map": "portal_gencode.gene.map",
  "input_checksum": "sha256:...",
  "input_genes": 10,
  "recognized_genes": "<computed>",
  "genes_with_coordinates": "<computed>",
  "genes_with_annotations": "<computed>",
  "unresolved_genes": ["PLINK1"],
  "input_gene_qc_file": "input_gene_qc.tsv",
  "validation_status": "PASS_WITH_WARNINGS",
  "execution_status": "COMPLETED",
  "exit_code": 0,
  "convergence_status": null,
  "num_chains": 10,
  "max_iterations": 500,
  "completed_iterations": 500,
  "final_iteration_index": 499,
  "log_file": "pigean_run.log",
  "timestamp": "2026-09-09T19:30:00Z"
}
```

> Fields use `null` when unavailable (e.g., `convergence_status` until Phase 5 implements trace analysis). The manifest distinguishes three statuses: `validation_status` (input QC), `execution_status` (engine exit), and `convergence_status` (model convergence).

This saves enormous headaches six months later when someone asks: *"Which gene-set database produced this result?"*

### Updated WDL Design

```wdl
workflow rock_pigean {
    input {
        # Required
        File input_file
        String output_name

        # User-facing options (not engine internals)
        String analysis_type = "positive-controls"
        String preset = "standard"
        String gene_set_profile = "default"
        String genome_build = "hg19"

        # Optional user files
        Array[File]? custom_gene_sets
        File? background_gene_list

        # Runtime
        Int memory_gb = 8
        Int cpu = 4
        Int disk_gb = 50

        # Developer/advanced only — not for normal users.
        # Allowing arbitrary image replacement undermines the golden test
        # and provenance guarantees. The manifest always records the actual digest.
        # String docker_image = "gcr.io/nitrogenase-docker/rock-pigean:2.0.0"
    }
    ...
}
```

Note: `gene_set_profile = "default"` resolves to the baked-in files inside the container. `custom_gene_sets` are actual WDL `File`s localized by the workflow engine. This avoids mixing container-internal paths with workflow inputs.

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| `run_pigean.py` core | 1–2 days | 1 day |
| Input validation | 0.5 day | 0.5 day |
| Run manifest | 0.5 day | Small |
| Basic text report (input QC, model status, ranked candidates, score definitions) | 0.5–1 day | 0.5 day |
| WDL parameterization | 0.5 day | 0.5 day |
| Golden test integration | — | Included above |
| **Total** | **~3.5–4 days** | **~2.5 days** |

---

## 7. Phase 2 — Custom Gene Sets + Background Genes

### Goal

Let users bring their own pathway databases and specify a background gene population.

### What Changes

1. **`--gene-sets custom`** with user-provided files
2. **`--background`** for specifying the eligible/background gene population
3. Gene set profile system: `default`, `msigdb-only`, `mouse-only`, `custom`

### Gene Set File Format

Tab-delimited, no header (engine format):
```
pathway_name    gene1    gene2    gene3    ...
```

> **GMT format support.** Standard GMT (as used by GSEA/MSigDB) has three fields: `name → description → genes`. The engine format is `pathway → genes`. If GMT files are accepted, the wrapper must explicitly convert and strip the description field — otherwise the description string enters the gene list as a spurious gene. Support standard GMT separately from the engine format, with explicit conversion.

Users can provide KEGG, Reactome, tissue-specific, disease-specific, or completely custom gene sets.

### Background Genes

> **The background is part of the scientific question, not a mere technical parameter.**
> Genes outside the positive list are **not** necessarily known biological negatives. Using "negative controls" wording is misleading.

Wires up `--positive-controls-all-in`, allowing the user to explicitly define the eligible/background gene population. The effect of background specification on enrichment calibration will be quantified during validation before recommendations are made (see §16).

**Three gene populations to distinguish explicitly:**

| Population | Definition |
|---|---|
| **Eligible to enter the input list** | Genes that could realistically have been selected (e.g., for a proteomics-derived list: proteins that passed QC and were eligible for selection) |
| **Represented in the model** | Genes present in the PIGEAN gene universe after reference/annotation intersection |
| **Eligible for the prediction report** | Genes that appear in ranked output |

Verify how `--positive-controls-all-in` affects each of these three populations rather than assuming it restricts all three.

Candidate background sets (appropriate choice depends on the analysis):
- For proteomics-derived input: **genes corresponding to proteins that passed QC and were eligible for selection** — then document mapping and exclusions
- Genes expressed in a tissue of interest
- A custom background matching the population of genes that could realistically have entered the analysis

> **"All protein-coding genes" should not silently become the universal default.** The background must match the experimental design.

### Interface

```bash
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --gene-sets custom \
    --custom-gene-set-files my_kegg.txt my_reactome.txt \
    --background all_protein_coding.txt \
    --output results/run2
```

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| Custom gene set wiring | Small | Small–medium |
| Background gene support | Small | Medium (verify enrichment calibration) |
| Gene set profiles | Small | Small |
| **Total** | **~1 day** | **~1–2 days** |

---

## 8. Phase 3 — Multi-Evidence Input Support

### Goal

Accept precomputed gene-level evidence — Bayes factors, Z-scores, percentiles, exome associations. This proves the evidence-adapter pattern before tackling the much harder GWAS adapter.

### Why Before GWAS

Conceptually, all evidence types ultimately contribute gene-level evidence to PIGEAN. The wrapper's responsibility is to validate and normalize each evidence type and invoke the corresponding `priors.py` interface; it does not necessarily materialize a common gene-level evidence file — `priors.py` may perform parts of that conversion internally.

Gene-level inputs (BF, Z, exome) are the simplest form of this evidence. They don't require coordinate mapping, genome builds, LD correction, or S2G assignment. Getting them working first means:

1. The adapter pattern is validated.
2. GWAS becomes "just another adapter" — the wrapper validates and normalizes the input, then delegates to the engine's own GWAS processing.

### Interface

```bash
# From precomputed Bayes factors — must declare scale explicitly
run_pigean.py \
    --analysis gene-scores \
    --input gene_bayes_factors.tsv \
    --input-format bayes-factors \
    --bf-scale log-natural \
    --output results/run3

# Supported --bf-scale values (v1.1):
#   raw          — raw Bayes factor
#   log-natural  — natural log BF (matches engine's expected log_bf column)
#   log10        — log10 BF
#
# Probability input is NOT a BF scale — it is a separate input format:
#   --input-format gene-probabilities
# The engine converts probabilities using:
#   log_bf = log(prob / (1 - prob)) - background_log_bf
# That conversion depends on background odds. An externally supplied
# posterior probability calculated with a different prior does not
# correspond to the same Bayes factor. Probability support requires
# documented prior assumptions matching the engine's conversion.

# From gene probabilities (separate format, not a BF scale)
run_pigean.py \
    --analysis gene-scores \
    --input gene_probabilities.tsv \
    --input-format gene-probabilities \
    --output results/run3b

# From exome associations
run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --output results/run4
```

> **`--bf-scale` is always required** when `--input-format bayes-factors` is specified. There is no silent default scale. This prevents the most common source of silently wrong results.

> **Probability is not simply another BF scale.** The engine converts probabilities using `log(p/(1-p)) - background_log_bf`, which depends on background odds. An externally supplied posterior probability may have been calculated using a different prior, so its numeric value alone does not establish the corresponding Bayes factor. For v1.1: initially support raw BF, log-natural BF, and log10 BF. Treat probability input as a separate, explicitly documented input type, and require its meaning and prior assumptions to match the supported conversion before accepting it.

> **Z-score support requires validation before advertising.** The engine's Z-score reader contains different interpretations, including treating values as raw log-odds and mapping them to probabilities. The selected conversion must be validated against known results before general Z-score support is offered. Z-score support ships in v1.1 after explicit validation, not in v1.

### Evidence Adapter Pattern

Each input type has an adapter that translates user data into what `priors.py` expects:

```
positive-controls adapter:  gene list → --positive-controls-in
gene-scores adapter:        BF/Z/pct → --gene-bfs-in / --gene-zs-in / --gene-percentiles-in
exome adapter:              gene+p+beta → --exomes-in + column flags
gwas adapter (Phase 4):     SNP-level → --gwas-in + build + S2G + correction flags
```

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| Gene BF adapter (with scale handling) | Small–medium | Medium (validate each scale conversion) |
| Gene Z-score adapter | Small–medium | **Medium–high** (validate reader interpretations against known results) |
| Gene percentile adapter | Small | Small |
| Exome adapter | Small–medium | Medium |
| Adapter pattern framework | Medium | — |
| **Total** | **~2–3 days** | **~3 days** |

---

## 9. Phase 4 — GWAS + Genome Build Reference Packs

### Goal

Accept GWAS summary statistics and support hg38 via versioned reference packs.

### Why This Is the Hardest Phase

GWAS introduces all of these concerns at once:

- SNP coordinates → need genome build
- Genome build → need matching loc files (gene, TSS, exon)
- Column detection (CHR, BP, P, BETA, SE, N, FREQ)
- Effect alleles and orientation
- LD and gene-size correction (`--correct-huge`)
- SNP-to-gene mapping (proximity, credible sets, custom S2G)
- Chromosome naming (chr1 vs 1, chrX vs 23)
- Missing values, duplicates, coordinate ranges
- Sample size handling

This is why we validate gene-level evidence first (Phase 3).

### Interface

```bash
run_pigean.py \
    --analysis gwas \
    --input gwas_summary_stats.tsv.gz \
    --genome-build hg38 \
    --output results/run5
```

Auto-detects columns. User can override:
```bash
    --gwas-chrom-col CHR \
    --gwas-pos-col BP \
    --gwas-p-col PVALUE \
    --gwas-beta-col EFFECT \
    --gwas-se-col STDERR \
    --gwas-n 50000
```

> **Build declaration is required, not inferred.** `--genome-build` is mandatory for GWAS analysis. Coordinate ranges can identify some invalid positions, but they cannot establish that a file is genuinely GRCh38. The wrapper validates reference compatibility and reports coordinate convention checks as confirmatory, not definitive.

### Reference Packs

Instead of loose files, create versioned reference packs:

```
references/
├── hg19/
│   ├── manifest.yaml
│   ├── gene.loc
│   ├── tss.loc
│   └── exons.loc
└── hg38/
    ├── manifest.yaml
    ├── gene.loc
    ├── tss.loc
    └── exons.loc
```

Each manifest:

```yaml
reference_name: PIGEAN-GRCh38-v1
genome_build: GRCh38
gene_annotation: GENCODE v44
created: 2026-09-10
gene_count: 20,096
coordinate_system: 1-based, inclusive
source: https://www.gencodegenes.org/human/release_44.html
```

And separately, annotation packs (symbol-based, annotation-version-dependent — see §13):

```
annotations/
├── msigdb_2026/
│   ├── manifest.yaml
│   └── msigdb_nohp.txt
└── mouse_phenotype_2024/
    ├── manifest.yaml
    └── mouse_2024.txt
```

Every run records:
```
Genome reference: PIGEAN-GRCh38-v1
Gene sets: MSigDB-2026 + MousePhenotype-2024
Gene map: gencode-v44
PIGEAN image: sha256:...
```

### Creating hg38 References

Only the three `.loc` files are genome-build-specific. Gene sets and gene maps use gene symbols (not build-specific, but annotation-version-dependent — see §13). For the current gene/TSS/exon mapping path, these three `.loc` files are the minimum GRCh38-specific assets. Additional GWAS/S2G workflows (Phase 4+) may require further build-matched reference resources (e.g., LD panels, S2G maps). hg38 support requires at minimum:

1. `GRCh38.gene.loc` — derived from GENCODE v44 GTF
2. `GRCh38.tss.loc` — TSS positions from same source
3. `GRCh38.exons.loc` — exon boundaries from same source

These can be created with a straightforward script against GENCODE annotations.

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| GWAS adapter + column detection | Medium | **High** |
| hg38 reference files | Medium | **High** (coordinate validation) |
| Reference pack system | Small–medium | Small |
| Build selection in wrapper | Small | Medium |
| Chromosome naming normalization | Small | Medium |
| **Total** | **~3–4 days** | **~3–5 days** |

---

## 10. Phase 5 — Full Interpretation Report and Convergence Diagnostics

> **Note:** A basic text report (input QC, model status, ranked candidates, score definitions) is delivered in Phase 1 as part of v1. Phase 5 adds the full HTML report, interactive visualizations, convergence diagnostics, and factor analysis summaries.

### Goal

Produce a full HTML report that makes the four raw output files interpretable without manual inspection. Extend the basic v1 text report with interactive visualizations, detailed convergence assessment, and pathway decomposition summaries.

### Interpretation Contract (Required Before Report Implementation)

> **Before generating any report, define every displayed field: its source column, scale, transformation, and biological meaning.**

PIGEAN's documentation distinguishes direct genetic support, annotation-derived support, and combined support; some displayed quantities are expressed on log-odds scales. Labels such as "prior" and "posterior" must be verified against the exact engine version before use in reports.

| Field Label | Source Column | Scale | Transformation | Biological Meaning |
|---|---|---|---|---|
| (To be populated during Phase 1 implementation — every field in the report must have an entry here before it is displayed) | | | | |

This table becomes a documented contract that the report code implements. No field is displayed without an entry.

### Why This Is Core, Not Optional

After manually interpreting the outputs during development, the answer is clearly yes — a report is essential. Scientists should not need to open `gs.out`, `gss.out`, `ggss.out`, and `p.out` separately and write custom `awk` commands to understand their results.

### Report Structure

> **Illustrative report layout.** Values marked with † are placeholders not extracted from a real run. Values from the actual 10-gene golden run (LEP prior=1.89, ADIPOQ=1.36, etc.) are used where available; all others are illustrative.
>
> **Important: input genes receive a default input probability of 0.95.** High output probabilities for these genes partly reflect this supplied assumption — they are not independent validation. The report must display input-positive status beside each gene to make this visible.

```
PIGEAN Analysis Report
══════════════════════════════════════════════════════

INPUT QC
────────
Analysis type:        positive-controls
Input genes:          10
Recognized symbols:   <value> (of 10 found in gene map)
With coordinates:     <value> (mapped to hg19 loci)
With ≥1 gene set:    <value> (have pathway annotations)
Unresolved genes:      1 (PLINK1) — see input_gene_qc.tsv
Gene sets retained:   2,094
Genome build:         hg19
Validation:           PASS_WITH_WARNINGS ⚠

MODEL STATUS
────────────
Engine status:        Completed ✓
Chains:               10
Iterations completed: 500 (index 0–499)
Convergence:          Not assessed
                      (reliable convergence assessment requires
                       trace analysis — see Phase 5)
Learned p:            0.00356
Learned σ²:           6.41e-08

  Note: "Completed" means the engine exited successfully.
  It does not imply convergence has been verified.

GENE PRIORITIZATION — Input Genes
──────────────────────────────────
Gene      Input?  Prior   Prior %ile †  Posterior   Top Pathway (from ggss.out)
LEP       yes     1.89    <value> †     0.986       mp_increased_pancreatic_islet_number (β=0.022)
ADIPOQ    yes     1.36    <value> †     0.979       <derived from ggss.out> †
BRD2      yes     0.80    <value> †     0.969       <derived from ggss.out> †
SOX2      yes     0.74    <value> †     0.967       <derived from ggss.out> †
LITAF     yes     0.22    <value> †     0.957       <derived from ggss.out> †
LECT2     yes     0.20    <value> †     0.955       <derived from ggss.out> †
SLCO1B1   yes     0.18    <value> †     0.955       <derived from ggss.out> †
RMI2      yes     0.05    <value> †     0.951       <derived from ggss.out> †
CD300LG   yes     0.04    <value> †     0.951       <derived from ggss.out> †
PLINK1    yes    -0.01    —             0.950       (unresolved — 0 gene-set memberships)

  Note: Input genes receive a default input probability of ≈0.95.
  High posteriors for input genes reflect this assumption and are
  not independent evidence of biological relevance.

ADDITIONAL PRIORITIZED CANDIDATES (not in input)
─────────────────────────────────────────────────
Gene      Input?  Prior   Annotation Coverage  Posterior
LEPR      no      0.988   677 gene sets        0.174
PPARG     no      0.640   1,003 gene sets      0.124
APOE      no      0.432   1,238 gene sets      0.097
IRS2      no      <value> † <value> †          <value> †
LDLR      no      <value> † <value> †          <value> †

WHY DOES LEPR SCORE HIGH?
──────────────────────────
(Associated gene sets and their joint effects from ggss.out)
Pathway                                         Beta (joint)  Membership
<populated from ggss.out at report time> †      <value> †     ✓
...

  Note: These are associated gene sets, not exact causal contributions.
  The model's scaling and adjustments affect displayed values.

TOP ENRICHED PATHWAYS
─────────────────────
Pathway                                         Beta (joint)  Beta_tilde (marginal)  P (marginal)  Source
mp_increased_pancreatic_islet_number            0.0224        <value> †              0.029         mouse_2024
mp_increased_interscapular_fat_pad_weight       0.0191        <value> †              0.028         mouse_2024
WP_LEPTIN_AND_ADIPONECTIN                       0.0094        <value> †              0.013         msigdb_nohp
...
```

### Convergence Diagnostics (Phase 5 — Full Assessment)

Use trace outputs (`--gene-set-stats-trace-out`, `--betas-trace-out`, `--gene-stats-trace-out`) internally to assess convergence. Report shows:

```
Convergence: PASS ✓        (or)        Convergence: WARNING ⚠
```

Scientists should not need to inspect trace files manually.

> **v1 reports "Convergence: Not assessed."** Full convergence assessment (trace analysis, R̂ computation) is a Phase 5 deliverable. A process finishing successfully does not automatically mean "Model QC: PASS." Where convergence cannot yet be reliably assessed, the report says "not assessed" — never "PASS."

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| HTML report generation | 2–3 days | 1 day |
| Convergence diagnostics | 1 day | 1 day |
| "Why does gene X score high?" logic | 0.5 day | 0.5 day |
| **Total** | **~4 days** | **~2–3 days** |

---

## 11. Phase 6 — Factor / PheWAS / Advanced Modes

### Goal

Expose pathway decomposition and phenome-wide association analysis for advanced users.

### What This Unlocks

The `factor` mode decomposes a gene list into distinct biological mechanisms (sub-pathways). This answers a different question than enrichment:

- **Enrichment (gibbs):** "Which pathways overlap my genes?"
- **Factoring:** "What are the distinct biological themes within my gene list?"

PheWAS projects results across many phenotypes simultaneously.

### Interface

```bash
# Factor analysis
run_pigean.py \
    --analysis factor \
    --input genes.txt \
    --output results/factor_run

# PheWAS
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --phewas-stats phewas_data.tsv \
    --output results/phewas_run
```

### Additional Outputs

- `factors.out` — factor loadings
- `factors_anchor.out` — anchor loadings
- `*_clusters.out` — gene/gene-set/pheno clusters
- `phewas_stats.out` — phenome-wide associations

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| Factor mode wiring | Medium | High (interpret outputs) |
| PheWAS wiring | Medium | High |
| Anchoring options | Medium | High |
| **Total** | **~3 days** | **~4 days** |

---

## 12. Phase 7 — Production Infrastructure

### Goal

Make the pipeline production-ready and reproducible.

### Changes

- **Pin Docker base image** to an immutable digest (`FROM ubuntu@sha256:...`), not a mutable release tag
- **Multi-arch Docker builds** (amd64 + arm64)
- **WDL `scatter`** for running multiple gene lists in parallel
- **Preemptible/spot instance** support in WDL runtime
- **Log capture** as workflow output
- **Cross-validation** via `--hold-out-chrom`
- **CI/CD:** golden test runs on every PR
- **Container registry migration:** GCR → Artifact Registry (GCR is deprecated)

**Estimated effort:**

| Sub-task | Implementation | Validation |
|---|---|---|
| Dockerfile hardening | Small | Small |
| CI/CD pipeline | 1 day | 0.5 day |
| Scatter support | 0.5 day | 0.5 day |
| Spot/preemptible | Small | Small |
| **Total** | **~2 days** | **~1–2 days** |

---

## 13. Reference Pack Strategy

### Current State (hg19 Only, Loose Files)

```
data/
├── gene_set_list_mouse_2024.txt       # Symbol-based (annotation-version-dependent; see note)
├── gene_set_list_msigdb_nohp.txt      # Symbol-based (annotation-version-dependent; see note)
├── portal_gencode.gene.map            # Symbol-based (annotation-version-dependent; see note)
├── NCBI37.3.plink.gene.loc            # hg19 ONLY
├── refGene_hg19_TSS.subset.loc        # hg19 ONLY
└── NCBI37.3.plink.gene.exons.loc      # hg19 ONLY
```

### Target State (Versioned Packs)

```
references/
├── hg19/
│   ├── manifest.yaml              # Version, source, gene count, build
│   ├── gene.loc
│   ├── tss.loc
│   └── exons.loc
├── hg38/
│   ├── manifest.yaml
│   ├── gene.loc
│   ├── tss.loc
│   └── exons.loc

annotations/
├── mouse_phenotype_2024/
│   ├── manifest.yaml              # Version, gene set count, source
│   └── gene_sets.txt
├── msigdb_nohp/
│   ├── manifest.yaml
│   └── gene_sets.txt

gene_maps/
├── portal_gencode.gene.map
└── gencode_v44.gene.map
```

### Key Insight

Only the three `.loc` files are genome-build-specific. Gene sets and gene maps use gene symbols, which are **not** genome-build-specific but **are** annotation-version-dependent — gene symbols can be added, retired, or aliased across GENCODE/NCBI releases. Pin the annotation version in each pack's manifest and version them together so a reference pack is a complete, reproducible unit.

---

## 14. Evidence Adapter Architecture

This is the central abstraction that makes generalization clean.

```
                    ┌──────────────────────┐
                    │   EVIDENCE LAYER     │
                    │                      │
    positive ───────┤                      │
    controls        │                      │
                    │     evidence         │
    GWAS ───────────┤     adapters         │──── Normalized / engine-ready
                    │                      │     evidence
    exome ──────────┤                      │
                    │                      │
    gene BF/Z/pct ──┤                      │
                    └──────────┬───────────┘
                               │
                    ┌──────────▼───────────┐
                    │  ANNOTATION LAYER    │
                    │                      │
                    │  MSigDB              │
                    │  Mouse phenotypes    │──── Gene set matrix X
                    │  Custom gene sets    │
                    │                      │
                    └──────────┬───────────┘
                               │
                    ┌──────────▼───────────┐
                    │  REFERENCE LAYER     │
                    │                      │
                    │  hg19 / hg38         │
                    │  Gene map            │──── Coordinates + ID mapping
                    │  (versioned packs)   │
                    │                      │
                    └──────────┬───────────┘
                               │
                    ┌──────────▼───────────┐
                    │  INFERENCE LAYER     │
                    │                      │
                    │  PIGEAN Gibbs        │──── priors.py (engine)
                    │  (priors.py)         │
                    │                      │
                    └──────────┬───────────┘
                               │
                    ┌──────────▼───────────┐
                    │  INTERPRETATION      │
                    │  LAYER               │
                    │                      │
                    │  gs / gss / ggss / p │
                    │  QC + manifest       │──── Scientist-ready output
                    │  HTML report         │
                    │  Convergence check   │
                    │                      │
                    └──────────────────────┘
```

Each layer is independently testable and swappable. The wrapper orchestrates them.

### Adapter Contract

Every evidence adapter implements the same interface. Critically, adapters return **normalized files plus mapping decisions, exclusions, and evidence semantics** — not only a path.

```python
@dataclass
class NormalizationResult:
    """What came out of normalization — the file and the audit trail."""
    normalized_path: str          # Path to engine-ready file
    evidence_scale: str           # e.g., "log_bf", "probability", "z_score", "binary"
    species: str                  # e.g., "human" — v1 supports human only
    identifier_type: str          # e.g., "gene_symbol", "ensembl_id"
    input_count: int              # Total input records
    mapped_count: int             # Records successfully mapped
    excluded: list[dict]          # [{gene: "PLINK1", reason: "not in gene map"}, ...]
    duplicates_resolved: list[dict]  # [{gene: "TP53", action: "kept first"}, ...]
    ambiguous_aliases: list[dict]    # [{input: "HER2", resolved_to: "ERBB2", status: "unambiguous"}, ...]
    custom_sets_action: str       # "replace" or "supplement" (for gene sets)
    warnings: list[str]           # Free-text warnings for manifest/report


class EvidenceAdapter:
    """Base contract for all evidence adapters."""

    def validate(self, input_path: str, params: dict) -> ValidationResult:
        """Check input file format, required columns, value ranges.
        Returns PASS / PASS_WITH_WARNINGS / FAIL + list of issues."""

    def normalize(self, input_path: str, params: dict) -> NormalizationResult:
        """Transform user-provided file into the format priors.py expects.
        Returns NormalizationResult with the normalized file AND all mapping decisions."""

    def build_engine_args(self, result: NormalizationResult, params: dict) -> list[str]:
        """Return the priors.py CLI arguments for this evidence type.
        Example: ['--positive-controls-in', '/tmp/normalized_genes.txt']"""

    def qc_summary(self, result: NormalizationResult, params: dict) -> dict:
        """Return a dict of QC metrics for the run manifest and report.
        Includes all mapping decisions, exclusions, and evidence semantics."""
```

> **v1 scope declaration:** Human genes only. Gene symbols as the primary identifier (with documented alias resolution). Explicit support for other species or identifier types is deferred. If a user provides non-human data, the wrapper should fail with a clear message rather than silently mis-mapping. **Human-only is a declared input requirement** — gene symbol lookup alone cannot reliably establish the species of every submitted gene.

> **Alias resolution policy:** Accept unambiguous mappings automatically (e.g., HER2 → ERBB2 when the mapping is 1:1). For ambiguous aliases (one input symbol maps to multiple current symbols), require the user to supply an explicit mapping file or exclude the gene with a documented reason. Recording a guess in the audit trail does not make it reliable.

> **Custom gene sets:** The adapter must explicitly declare whether user-provided sets **replace** or **supplement** the defaults. This decision is surfaced in the run manifest and report.

This means adding a new evidence type (e.g., credible sets in Phase 4+) requires only a new adapter class — no changes to the wrapper's orchestration logic.

---

## 15. User Interface Design — Normal vs. Advanced

### Normal Interface (What Scientists See)

```bash
# What kind of analysis?
# What's your data?
# Which pathways to use?
# Optional: background genes
# Which genome build?
# How thorough?
# Where to put results?
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --gene-sets default \
    --background protein_coding.txt \
    --genome-build hg19 \
    --preset standard \
    --output results/my_run
```

That's it. Seven parameters, all with sensible defaults except `--input`.

### Advanced Interface (Power Users)

```bash
run_pigean.py \
    --analysis positive-controls \
    --input genes.txt \
    --config advanced.yaml
```

Where `advanced.yaml`:
```yaml
num_chains: 15
max_num_iter: 1000
prune_gene_sets: 0.7
sparse_frac: 0.02
convergence_threshold: 1.005
max_frac_sem: 0.005
hold_out_chrom: 6
```

Every run produces `resolved_config.yaml` — the full merged configuration (defaults + overrides) that was actually used. This is the single source of truth for reproducing the run:

```yaml
# resolved_config.yaml — auto-generated, do not edit
analysis: positive-controls
preset: standard
genome_build: hg19
num_chains: 15          # overridden from advanced.yaml (default: 10)
max_num_iter: 1000      # overridden from advanced.yaml (default: 500)
prune_gene_sets: 0.7    # overridden from advanced.yaml (default: 0.8)
sparse_frac: 0.02       # overridden from advanced.yaml (default: 0.01)
gene_set_filter_value: 0.01
gene_filter_value: 1
debug_level: 3
...
```

These parameters are documented but not surfaced in the standard help text. Users should not need to make decisions about chain counts or sparsity fractions unless they have a specific reason.

### What NOT to Put in the Normal Interface

| Parameter | Why Not |
|---|---|
| `num_chains` | 10 is the current pipeline default; keep hidden unless a use case requires override |
| `max_num_iter` | Convergence-based stopping handles this |
| `sparse_frac` | Engine internals |
| `sigma2` | Learned automatically |
| `prune_gene_sets` | 0.8 is the current engine default; alternative values remain advanced until benchmarked |
| `debug_level` | Users want results, not debug logs |
| `batch_size` | Memory optimization detail |
| `max_frac_sem` | Convergence detail |

---

## 16. Known Warnings — Verified vs. Hypothesized

> **Rule:** If we have not verified a behavior in code or benchmarked it, label it as a hypothesis — not a fact.

### Verified ✓

| Warning | Status | Evidence |
|---|---|---|
| `PLINK1` is unresolved | **Verified: gene has missing coordinates and zero annotation memberships.** N=0, Chrom=NA, coordinates missing, 0 gene-set memberships. The wrapper should report which individual lookup step failed (gene-map recognition, coordinate mapping, gene-set membership). The exact failed lookup step will be confirmed during Phase 0 implementation. | Observed in output: `PLINK1 N=0 Chrom=NA Start=NA End=NA` |

### Needs Verification ⚠ (Reclassified)

| Warning | Status | What We Need to Do |
|---|---|---|
| `--correct-huge was not used` | **Reclassified: needs verification.** The warning occurs in `_correct_beta_tildes()` when correction information is missing. The parser exposes `--no-correct-huge` (correction enabled by default). Original assessment that this is "harmless in gene-list mode" requires verification that the code path is truly inert without GWAS input. | Trace the `_correct_beta_tildes()` call path with gene-list-only input. Verify whether the correction logic is skipped entirely or executes with missing data. |
| `random.seed(0)` establishes reproducibility | **Hypothesis: reproducibility is incomplete.** The engine calls `random.seed(0)`, but also uses NumPy random sampling without a corresponding `np.random.seed()` in the inspected file. This does **not** establish reproducible Gibbs sampling. | Run the golden test case 5+ times and measure output variability. If variability is observed, the golden test must use statistical tolerances, not exact matching. |

### Hypothesized — Needs Verification ⚠

| Warning | Hypothesis | What We Need to Do |
|---|---|---|
| `Specified positive controls without --positive-controls-all-in` | **Hypothesis: this inflates enrichment scores.** The warning says it "may result in inflated enrichments" but we have not quantified the effect or verified the mechanism. | Run the same 10-gene list with and without `--positive-controls-all-in` (using all protein-coding genes as background). Compare priors and pathway betas. |
| `A large fraction of betas (0.0175) are of opposite signs` | **Hypothesis: fractions < 0.02 are "normal."** We have not benchmarked what fraction indicates a problem. | Run multiple gene lists of varying sizes. Record the opposite-sign fraction. Establish empirical thresholds. |
| Pruning recommendation (0.8 → 0.7 → 0.6) | **Hypothesis: lower pruning helps with collinear gene sets.** Plausible but not tested. | Benchmark with different `--prune-gene-sets` values on a collinear gene set case. |
| `naive_priors` is a good "fast" preset | **Hypothesis: naive_priors gives similar gene rankings to gibbs.** We have not compared them. | Run both on the golden test case. Compare rank correlation of gene priors. |

These hypotheses should be tested during Phase 0 or Phase 1 before they appear in user documentation.

---

## 17. Effort Estimates — Implementation vs. Validation

> "Wire a flag" ≠ "production-quality support."

| Phase | Implementation | Validation | Total |
|---|---|---|---|
| **Phase 0:** Golden test + repeated-run baselines | 2–3 hours | 0.5 day (repeated runs) | **~1 day** |
| **Phase 1:** Wrapper + WDL + validation + manifest + **basic report** | 3–4 days | 2 days | **~6 days** |
| **Phase 2:** Custom gene sets + background | 1 day | 1–2 days | **~2–3 days** |
| **Phase 3:** Multi-evidence input adapters | 2–3 days | 3 days | **~5–6 days** |
| **Phase 4:** GWAS + hg38 reference packs | 3–4 days | 3–5 days | **~7–9 days** |
| **Phase 5:** Full HTML report + convergence diagnostics | 4 days | 2–3 days | **~6–7 days** |
| **Phase 6:** Factor / PheWAS | 3 days | 4 days | **~7 days** |
| **Phase 7:** Production infra | 2 days | 1–2 days | **~3–4 days** |
| | | | |
| **Total** | **~20 days** | **~17–22 days** | **~37–43 days** |

> **Timeline note:** These totals sum to approximately **37–43 working days**, or **7.5–9 working weeks** for one person, before contingency. The **two-week v1 prototype target** (Phases 0–2, ~9 days) is realistic. The full schedule is provisional — treat per-phase estimates as independent commitments, not a single deadline.

### Key Insight on Validation Costs

| Feature | Wire it | Validate it |
|---|---|---|
| Gene BF/Z input | Small | Small–medium |
| Custom gene sets | Small | Medium |
| Exome adapter | Small–medium | Medium |
| GWAS adapter | Medium | **High** |
| GRCh38 references | Medium | **High** |
| Factor mode | Medium | High |

GWAS and hg38 have disproportionately high validation costs because getting coordinates wrong means scientifically incorrect results, silently.

---

## 18. MVP Definitions

### Revised Release Order

| Release | Scope | Phases |
|---|---|---|
| **v1: Usable gene-list workflow** | Phases 0–2, basic text report, pinned existing reference bundle, logs, automated regression checks, installation guide, and example dataset | 0, 1, 2 + parts of 7 |
| **v1.1: Validated gene evidence** | Explicit BF/log-BF support (validated), followed by separately validated Z-score and exome inputs | 3 (incremental) |
| **v2: GWAS support** | Validated GWAS processing, GRCh38 references, and compatibility checks | 4 |
| **v2.1: Full HTML report + convergence** | Full HTML interpretation report, convergence diagnostics, interactive visualizations | 5 |
| **Later releases** | Factor/PheWAS modes, batch execution, multi-arch builds, and additional platforms | 6, 7 |

### v1 — Usable Gene-List Workflow (Phases 0–2 + Basic Report)

A user can:
1. Provide a positive-control gene list ✓
2. Use default or custom gene sets ✓
3. Specify an appropriate background gene population ✓
4. Pipeline validates gene IDs and annotation coverage ✓
5. Default settings reproduce the golden test ✓
6. Pipeline produces gs/gss/ggss/p ✓
7. Pipeline produces a **basic text report** with: input QC, model status (with "not assessed" for unverified convergence), ranked candidates (with input-positive status), and score definitions ✓
8. Pipeline produces a QC summary and run manifest ✓
9. Every run records **pinned engine version, Python/dependency versions, reference checksums, and container digest** ✓
10. **Captured logs, exit status, and automated regression checks** ✓
11. Separate statuses for input validation, execution completion, and convergence ✓
12. **Installation guide and example dataset** included ✓

**v1 acceptance criterion:**

> **Another researcher can install it, run their own supported gene list, understand every exclusion and displayed score, and reproduce the analysis without your help.**
>
> Have one researcher try that before expanding to additional modes.

**Timeline:** ~2 weeks (reasonable prototype target, provided the existing run is reproducible and references are ready)

### v1.1 — Validated Gene Evidence (Phase 3, incremental)

Everything in v1, plus:
1. Accept precomputed Bayes factors with explicit scale declaration (raw BF, log-natural BF, log10 BF) ✓
2. Accept gene probabilities as a separate input format (with documented prior assumptions) ✓
3. Accept gene Z-scores (after validating the engine's reader interpretations against known results) ✓
4. Accept exome associations ✓

> BF support (raw, log-natural, log10) ships first because the engine's expected format is well-defined. Probability, Z-score, and exome support ship only after their evidence conversions have been independently validated.

### v2 — GWAS Support (Phase 4)

Everything in v1.1, plus:
1. Accept GWAS summary statistics ✓
2. Support hg38 via reference packs ✓
3. GWAS column auto-detection and validation ✓

**Timeline:** ~4 weeks after v1

### v2.1 — Full HTML Report + Convergence (Phase 5)

Everything in v2, plus:
1. Full HTML interpretation report ✓
2. Automatic convergence diagnostics (trace analysis, R̂) ✓
3. Interactive visualizations ✓

### Later Releases (Phases 6–7)

Everything in v2, plus:
1. Factor / pathway decomposition mode ✓
2. PheWAS mode ✓
3. CI/CD with golden tests ✓
4. Scatter for parallel runs ✓
5. Spot/preemptible support ✓
6. Multi-architecture Docker builds ✓

**Timeline:** ~2 weeks after v2

> **Note:** The two-week v1 estimate is a reasonable prototype target. The full schedule totals approximately 37–43 working days (7.5–9 weeks) for one person, before contingency. The schedule is provisional until the evidence conversions and GWAS reference work have been validated.

---

## 19. Open Questions

### Decisions Needed Before Phase 1

1. **Wrapper language:** Python (matches priors.py) or bash (simpler)? Recommendation: Python — it can parse manifests, validate gene names against the gene map, generate JSON, and generate the HTML report.

2. **Docker registry:** Stay on GCR (`gcr.io/nitrogenase-docker/`) or migrate to Artifact Registry? GCR is deprecated.

3. **Container tagging:** ~~Should the wrapper version and priors.py version be tracked separately?~~ **Resolved:** Yes — the manifest now tracks `wrapper_version` and `engine.version` (git SHA or file checksum) separately, reflecting the engine-freeze architecture.

4. **Interpretation field definitions:** Before implementing any report output, define every displayed field (source column, scale, transformation, biological meaning) and verify labels like "prior" and "posterior" against the exact engine version. See §10 Interpretation Contract.

5. **Species and identifier scope:** v1 declares human support only with gene symbols as the primary identifier. Document alias resolution rules (e.g., HER2 → ERBB2). Should the wrapper hard-fail on non-human identifiers?

### Decisions Needed Before Phase 2

6. **GMT format support:** ~~Should the wrapper accept standard GMT files?~~ **Decided: Yes.** Accept standard GMT files (`name → description → genes`) in addition to the engine format (`pathway → genes`). Implement explicit conversion that strips the description field to prevent it from entering the gene list.

7. **Custom gene set interaction:** Do user-provided gene sets **replace** or **supplement** the defaults? This must be explicitly declared and surfaced in the manifest.

### Decisions Needed Before Phase 3

8. **Bayes factor scale validation:** The engine expects `Gene` and `log_bf` (natural-log BF) columns. Validate the conversion path for each supported scale (raw BF, log10 BF) against known results before shipping.

9. **Probability input conversion:** The engine converts probabilities using `log(p/(1-p)) - background_log_bf`, which depends on background odds. Document required prior assumptions and validate that the conversion produces correct results before accepting probability inputs. This is a separate input format, not a BF scale.

10. **Z-score reader interpretation:** The engine's Z-score reader contains multiple interpretations. Which is correct for each use case? Validate before advertising support.

### Decisions Needed Before Phase 4

11. **hg38 gene annotation source:** GENCODE v44? NCBI RefSeq? Which version to standardize on?

12. **GWAS column auto-detection:** How aggressive should it be? priors.py already has some column inference. Should the wrapper add more?

13. **GWAS build declaration:** Require a declared build rather than inferring from coordinate ranges. Coordinate ranges can identify some invalid positions but cannot establish that a file is genuinely GRCh38.

### Decisions Needed Before Phase 6

14. **Factor mode as separate workflow or mode switch?** Factor analysis has different inputs (anchoring), outputs (factors, clusters), and parameters (phi, alpha0). May be cleaner as a separate WDL workflow.

15. **Which additional gene set databases to ship?** Candidates: KEGG, Reactome, GO (standalone), DisGeNET, GWAS Catalog, GTEx tissue-specific.

### Standing Questions

16. **Background gene validation:** Need to experimentally verify the impact of `--positive-controls-all-in` on enrichment calibration before documenting it as guidance. Verify how `--positive-controls-all-in` affects the three gene populations (eligible-to-enter, represented-in-model, eligible-for-report).

17. **Convergence thresholds:** What convergence R metric and SEM values actually indicate a problem? Verify whether the engine's reported quantity is specifically Gelman–Rubin R̂ before labeling it as such. Need empirical benchmarks, not assumptions.

18. **Fast preset validation:** Does `naive_priors` produce similar rankings to `gibbs`? Need head-to-head comparison on the golden test case before advertising it as a "fast" option.

19. **Input gene default probability:** Input genes receive a default probability of ≈0.95. Verify the exact mechanism and ensure the report clearly distinguishes model-supported findings from supplied assumptions. An unresolved input gene retaining ≈0.95 must not look like a confidently supported biological result.

---

## How to Run — Manual Commands for Each Phase

> All commands assume you are in the repo root: `cd /fsx/home/l112724/gwas_pipeline/rock-pigean`

### Prerequisites

```bash
# Python 3.x with numpy and scipy
python3 --version    # 3.11+ recommended
pip install numpy scipy   # if not already installed
```

### Source File Map — Which Code Lives Where

```
rock-pigean/                          ← repo root
│
├── priors.py                         [16,567 lines] PIGEAN engine — NEVER MODIFIED
│                                     The Gibbs sampler. 309 optparse flags. Treated as
│                                     a frozen dependency. SHA-256 locked in Phase 0.
│
├── run_pigean.py                     [372 lines] CLI entry point (Phase 1+2)
│                                     Parses --input/--output/--gene-sets/--background etc.,
│                                     orchestrates all pigean/ modules, exits 0 or 1.
│
├── pigean/                           Python package — the wrapper layer
│   ├── __init__.py                   [3 lines]   Version string ("1.0.0")
│   ├── config.py                     [106 lines] DEFAULTS dict, resolve_config() precedence
│   │                                              (DEFAULTS → JSON config → CLI args)
│   ├── references.py                 [184 lines] Resolve data file paths (Docker vs local),
│   │                                              gene-set profiles (default/mouse-only/msigdb-only/custom),
│   │                                              resolve_gene_set_paths() with custom+action support
│   ├── validation.py                 [307 lines] Gene QC: load_gene_map(), validate_input_genes(),
│   │                                              validate_background(), compute_validation_status(),
│   │                                              write_gene_qc_tsv()
│   ├── engine.py                     [153 lines] build_priors_command() → exact flag list,
│   │                                              save_command(), run_engine() via subprocess
│   ├── manifest.py                   [244 lines] generate_manifest() → run_manifest.json with
│   │                                              checksums, timing, three-status model, provenance
│   ├── report.py                     [184 lines] generate_report() → human-readable report.txt
│   │                                              (QC, execution, output stats, interpretation notes)
│   └── adapters/
│       ├── __init__.py               [92 lines]  ValidationStatus enum, dataclasses:
│       │                                          FileValidationResult, NormalizationResult,
│       │                                          GeneQCResult, GeneSetPreparationResult,
│       │                                          BackgroundQCResult, EvidenceAdapter base class,
│       │                                          get_adapter() registry, SUPPORTED_ANALYSIS_TYPES
│       ├── base.py                   [~130 lines] BaseEvidenceAdapter abstract class,
│       │                                          EvidenceNormalizationResult dataclass,
│       │                                          validate_tabular_file() shared utility
│       ├── positive_controls.py      [182 lines] PositiveControlsAdapter: validate_file(),
│       │                                          normalize() (input/ + normalized/), build_engine_args()
│       │                                          Also used for background gene normalization.
│       ├── bayes_factor.py           [~140 lines] BayesFactorAdapter: Gene/log_bf validation,
│       │                                          NA removal, duplicate warnings, --gene-bfs-in
│       ├── zscore.py                 [~110 lines] ZScoreAdapter: requires gene_column + score_column,
│       │                                          log-odds semantics, --gene-zs-in/id-col/value-col
│       ├── percentile.py            [~110 lines] PercentileAdapter: requires gene_column + score_column,
│       │                                          higher_is_better flag, --gene-percentiles-in
│       ├── exome.py                  [~120 lines] ExomeAdapter: minimal validation (engine auto-detects),
│       │                                          optional column overrides, --exomes-in
│       ├── gwas.py                   [~210 lines] GwasAdapter: GWAS summary stats validation,
│       │                                          column override mapping, --gwas-in + 11 override flags
│       │                                          Minimal normalization (no data transformation)
│       └── gene_sets.py              [275 lines] validate_gene_set_file(), convert_gmt_to_engine_format(),
│                                                  prepare_custom_gene_sets() — legal combinations gatekeeper
│
├── data/                             Reference data files (baked into Docker image)
│   ├── gene_set_list_mouse_2024.txt  Mouse phenotype gene sets (annotation DB)
│   ├── gene_set_list_msigdb_nohp.txt MSigDB gene sets (annotation DB)
│   ├── portal_gencode.gene.map       Gene symbol → ID mapping
│   ├── NCBI37.3.plink.gene.loc       hg19 gene coordinates
│   ├── refGene_hg19_TSS.subset.loc   hg19 TSS coordinates
│   └── NCBI37.3.plink.gene.exons.loc hg19 exon coordinates
│
├── ex/
│   └── gene_list                     Example 10-gene input (same as golden test input)
│
├── tests/
│   ├── test_config.py                [114 lines] 10 tests: defaults, precedence, Phase 2 config keys
│   ├── test_validation.py            [242 lines] 21 tests: file validation, normalization, gene QC,
│   │                                              combined status (Phase 1 + Phase 2)
│   ├── test_engine_command.py        [255 lines] 14 tests: command structure, flag values, background flag,
│   │                                              evidence adapter engine args
│   ├── test_gene_sets.py             [370 lines] 36 tests: GMT conversion, file validation,
│   │                                              legal combinations, resolve paths, background validation
│   ├── test_evidence_adapters.py     [271 lines] 47 tests: adapter registry, BF/zscore/percentile/exome
│   │                                              adapters (validate, normalize, engine args), backward compat
│   ├── test_gwas_adapter.py          [~330 lines] 45 tests: GwasAdapter (validate, normalize, engine args),
│   │                                              genome build normalization, hg38 reference failure,
│   │                                              GWAS engine command, config keys
│   ├── fixtures/
│   │   ├── evidence/                 6 fixture files: bayes_factor_valid.tsv, bayes_factor_invalid_nocol.tsv,
│   │   │                              bayes_factor_with_na.tsv, zscore_valid.tsv, percentile_valid.tsv, exome_valid.tsv
│   │   └── gwas/                     5 fixture files: valid_gwas_hg19.tsv, gwas_custom_cols.tsv,
│   │                                  gwas_missing_cols.tsv, gwas_no_header.txt, gwas_with_blanks.tsv
│   └── golden/
│       ├── run_golden_test.sh        [383 lines] Regression test runner (56 tolerance checks)
│       ├── checkpoints.json          LOCKED — stochastic tolerances from 6 baseline runs
│       ├── input/
│       │   └── positive_controls.txt Frozen 10-gene test input
│       └── expected/
│           ├── gs.golden.out         Golden gene stats output
│           ├── gss.golden.out        Golden gene-set stats output
│           ├── ggss.golden.out       Golden gene-gene-set stats output
│           └── p.golden.out          Golden model parameters output
│
├── wdl/
│   └── rock_pigean.wdl               [87 lines] WDL workflow: parameterized inputs including
│                                                  custom gene sets + background (Phase 2)
│
├── Dockerfile                         [17 lines] ADD data/, priors.py, run_pigean.py, pigean/
│                                                  ENTRYPOINT → run_pigean.py
│
├── docs/
│   ├── evidence_inputs.md             User-facing documentation for all 5 gene-level evidence types
│   └── gwas_input.md                  User-facing documentation for GWAS input + genome builds
│
├── plan.README.md                     This file — architecture plan + changelogs
└── README.md                          User-facing usage docs
```

#### Which files belong to which phase

| Phase | Production code | Tests | Infrastructure |
|-------|----------------|-------|----------------|
| **Phase 0** | `priors.py` (frozen) | `tests/golden/run_golden_test.sh`, `tests/golden/checkpoints.json`, `tests/golden/input/`, `tests/golden/expected/` | — |
| **Phase 1** | `run_pigean.py`, `pigean/__init__.py`, `pigean/config.py`, `pigean/references.py`, `pigean/validation.py`, `pigean/engine.py`, `pigean/manifest.py`, `pigean/report.py`, `pigean/adapters/__init__.py`, `pigean/adapters/positive_controls.py` | `tests/test_config.py`, `tests/test_validation.py`, `tests/test_engine_command.py` | `Dockerfile`, `wdl/rock_pigean.wdl`, `README.md` |
| **Phase 2** | `pigean/adapters/gene_sets.py` (new), modifications to `run_pigean.py`, `pigean/adapters/__init__.py`, `pigean/adapters/positive_controls.py`, `pigean/references.py`, `pigean/config.py`, `pigean/validation.py`, `pigean/engine.py`, `pigean/manifest.py`, `pigean/report.py` | `tests/test_gene_sets.py` (new), additions to `tests/test_config.py`, `tests/test_engine_command.py`, `tests/test_validation.py` | Updates to `wdl/rock_pigean.wdl` |
| **Phase 3** | `pigean/adapters/base.py` (new), `pigean/adapters/bayes_factor.py` (new), `pigean/adapters/zscore.py` (new), `pigean/adapters/percentile.py` (new), `pigean/adapters/exome.py` (new), modifications to `run_pigean.py`, `pigean/adapters/__init__.py`, `pigean/adapters/positive_controls.py`, `pigean/config.py`, `pigean/engine.py`, `pigean/validation.py`, `pigean/manifest.py`, `pigean/report.py` | `tests/test_evidence_adapters.py` (new), `tests/fixtures/evidence/` (new), rewritten `tests/test_engine_command.py` | Updates to `wdl/rock_pigean.wdl`, `Dockerfile`, `docs/evidence_inputs.md` (new) |
| **Phase 4** | `pigean/adapters/gwas.py` (new), modifications to `run_pigean.py`, `pigean/adapters/__init__.py`, `pigean/references.py`, `pigean/config.py`, `pigean/manifest.py`, `pigean/report.py` | `tests/test_gwas_adapter.py` (new), `tests/fixtures/gwas/` (new) | Updates to `wdl/rock_pigean.wdl`, `Dockerfile`, `docs/gwas_input.md` (new) |

#### Data flow: which file does what at runtime

```
run_pigean.py                         ← CLI parsing, orchestration
    │
    ├── pigean/config.py              ← Merge DEFAULTS + JSON + CLI → resolved config
    ├── pigean/references.py          ← Find data/ files for the chosen profile + build
    │
    ├── pigean/adapters/
    │   ├── positive_controls.py      ← Validate & normalize input gene file
    │   │                                (also normalizes background file if --background)
    │   └── gene_sets.py              ← Validate & convert custom gene-set files
    │                                    (GMT→engine format, legal combination checks)
    │
    ├── pigean/validation.py          ← Gene QC against gene map + loc files + gene sets
    │                                    Background QC (input genes in background?)
    │                                    Combined validation status
    │
    ├── pigean/engine.py              ← Build priors.py command → run subprocess
    │       │
    │       └── priors.py             ← Gibbs sampler (engine) → gs/gss/ggss/p.out
    │
    ├── pigean/manifest.py            ← Write run_manifest.json (full provenance)
    └── pigean/report.py              ← Write report.txt (human-readable summary)
```

---

### Phase 0 — Run the Golden Regression Test

**Purpose:** Verify the PIGEAN engine produces results within the locked stochastic tolerances (56 checks).

**Input required:**
- `tests/golden/input/positive_controls.txt` (10-gene list, already in repo)
- All data files in `data/` (already in repo)

**Option A — Run engine fresh and validate:**

```bash
bash tests/golden/run_golden_test.sh
```

This runs `priors.py` with the frozen input, then validates the output against `tests/golden/checkpoints.json`. Output goes to a temp directory and is cleaned up after.

```bash
# Keep the output directory for inspection:
bash tests/golden/run_golden_test.sh --keep
```

**Option B — Validate existing outputs:**

```bash
# Point at any directory containing gs.out, gss.out, ggss.out, p.out:
bash tests/golden/run_golden_test.sh /path/to/existing/output
```

**Expected result:** `56 PASS, 0 FAIL`

---

### Phase 1 — Run the Wrapper (Default Settings)

**Purpose:** Run the full pipeline through the wrapper layer with default parameters. Reproduces the golden test via `run_pigean.py`.

**Input required:**
- A gene list file (one gene symbol per line)

**Minimal run (default settings — matches Phase 0):**

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --output runs/my_run
```

**Explicit settings (same as above, all defaults spelled out):**

```bash
python3 run_pigean.py \
    --analysis positive-controls \
    --input ex/gene_list \
    --gene-sets default \
    --genome-build hg19 \
    --preset standard \
    --output runs/my_run_explicit
```

**With a JSON config override:**

```bash
# Create a config file:
echo '{"max_num_gene_sets": 3000, "debug_level": 2}' > my_config.json

python3 run_pigean.py \
    --input ex/gene_list \
    --config my_config.json \
    --output runs/my_run_custom_config
```

**Overwrite an existing output directory:**

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --output runs/my_run \
    --overwrite
```

**Output structure:**

```
runs/my_run/
├── input/
│   └── positive_controls.txt       # Copy of original input
├── normalized/
│   └── positive_controls.txt       # Normalized (deduped, stripped)
├── gs.out                          # Gene stats (42,216 rows)
├── gss.out                         # Gene-set stats (~41,627 rows)
├── ggss.out                        # Gene-gene-set stats
├── p.out                           # Model parameters
├── input_gene_qc.tsv              # Per-gene QC table
├── priors_command.txt             # Exact engine command used
├── resolved_config.json           # Final merged config
├── run_manifest.json              # Full provenance
├── pigean_run.log                 # Engine stdout/stderr
└── report.txt                     # Human-readable summary
```

**Validate the wrapper output against golden test:**

```bash
bash tests/golden/run_golden_test.sh runs/my_run
```

---

### Phase 2 — Custom Gene Sets

**Purpose:** Run with built-in gene-set profiles or user-supplied gene-set files.

**Input required:**
- A gene list file
- (Optional) Custom gene-set file(s) in engine or GMT format

#### Built-in profiles

**Mouse phenotype gene sets only (1 database instead of 2):**

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets mouse-only \
    --output runs/mouse_only
```

**MSigDB gene sets only:**

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets msigdb-only \
    --output runs/msigdb_only
```

#### Custom gene sets — engine format

Engine format is tab-delimited, no header: `pathway_name<TAB>gene1<TAB>gene2<TAB>...`

```bash
# Create a sample custom gene-set file:
cat > /tmp/my_pathways.txt << 'EOF'
LIPID_METABOLISM	LEP	ADIPOQ	PPARG	LEPR	PCSK9
IMMUNE_SIGNALING	BRD2	LITAF	CD300LG	SOX2
EOF

# Replace all built-in gene sets with your custom file:
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets custom \
    --custom-gene-set-files /tmp/my_pathways.txt \
    --custom-gene-set-format engine \
    --custom-gene-set-action replace \
    --output runs/custom_replace

# Supplement the default gene sets with your custom file:
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets default \
    --custom-gene-set-files /tmp/my_pathways.txt \
    --custom-gene-set-format engine \
    --custom-gene-set-action supplement \
    --output runs/custom_supplement
```

#### Custom gene sets — GMT format

GMT format (as used by GSEA/MSigDB): `pathway_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...`

The wrapper automatically strips the description column during conversion.

```bash
# Create a sample GMT file:
cat > /tmp/my_pathways.gmt << 'EOF'
LIPID_METABOLISM	Lipid metabolism pathway	LEP	ADIPOQ	PPARG	LEPR
IMMUNE_SIGNALING	Immune signaling genes	BRD2	LITAF	CD300LG
EOF

# Use GMT with supplement (adds to default gene sets):
python3 run_pigean.py \
    --input ex/gene_list \
    --custom-gene-set-files /tmp/my_pathways.gmt \
    --custom-gene-set-format gmt \
    --custom-gene-set-action supplement \
    --output runs/gmt_supplement
```

#### Multiple custom gene-set files

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets custom \
    --custom-gene-set-files /tmp/pathways_a.txt /tmp/pathways_b.txt \
    --custom-gene-set-format engine \
    --custom-gene-set-action replace \
    --output runs/multi_custom
```

---

### Phase 2 — Background Genes

**Purpose:** Specify the eligible gene population (background) for enrichment analysis.

**Input required:**
- A gene list file (input genes)
- A background gene list file (one gene per line — should be a superset of the input genes)

```bash
# Create a sample background file (superset of input genes):
cat > /tmp/background_genes.txt << 'EOF'
LEP
ADIPOQ
BRD2
SOX2
LITAF
SLCO1B1
RMI2
LECT2
CD300LG
LEPR
PPARG
APOE
PCSK9
EOF

# Run with background:
python3 run_pigean.py \
    --input ex/gene_list \
    --background /tmp/background_genes.txt \
    --output runs/with_background
```

The wrapper maps `--background` to the engine's `--positive-controls-all-in` flag. The report and manifest will show:
- How many background genes were recognized
- Which input genes are missing from the background (warning)

#### Combined: custom gene sets + background

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets mouse-only \
    --custom-gene-set-files /tmp/my_pathways.gmt \
    --custom-gene-set-format gmt \
    --custom-gene-set-action supplement \
    --background /tmp/background_genes.txt \
    --output runs/full_phase2
```

---

### Phase 3 — Multi-Evidence Input Types

**Purpose:** Run PIGEAN with gene-level evidence beyond positive controls — Bayes factors, Z-scores, percentiles, or exome associations.

**Input required:**
- A tab-separated evidence file with gene identifiers and scores/statistics

#### Gene-Level Bayes Factors

The simplest evidence adapter. Input file must have `Gene` and `log_bf` columns (natural-log Bayes factors):

```bash
# Using default column names (Gene, log_bf):
python3 run_pigean.py \
    --analysis gene-bayes-factor \
    --input evidence.tsv \
    --output runs/bf_run
```

#### Gene-Level Z-Scores (Log-Odds Mode)

**Requires** `--gene-column` and `--score-column`. The engine treats these as log-odds, NOT statistical Z-scores:

```bash
python3 run_pigean.py \
    --analysis gene-z-score \
    --input scores.tsv \
    --gene-column Gene \
    --score-column Score \
    --output runs/zscore_run
```

> **Known engine bug:** Z-score mode may produce a scipy domain error with certain input distributions. This is a frozen-engine issue, not a wrapper bug.

#### Gene-Level Percentiles

**Requires** `--gene-column` and `--score-column`:

```bash
# Lower score = stronger evidence (default)
python3 run_pigean.py \
    --analysis gene-percentile \
    --input ranks.tsv \
    --gene-column Gene \
    --score-column rank_score \
    --output runs/percentile_run

# Higher score = stronger evidence
python3 run_pigean.py \
    --analysis gene-percentile \
    --input enrichment.tsv \
    --gene-column Gene \
    --score-column enrichment_score \
    --higher-is-better \
    --output runs/percentile_higher
```

> **Known engine bug:** Percentile mode crashes on `options.top_posterior` attribute lookup (engine line 16343). Frozen-engine issue.

#### Exome Associations

Engine auto-detects columns. Requires ≥2 of {p-value, beta, SE}:

```bash
# Minimal (auto-detect columns):
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --output runs/exome_run

# With column overrides:
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --exomes-gene-col GENE_SYMBOL \
    --exomes-p-col PVALUE \
    --exomes-n 50000 \
    --output runs/exome_custom
```

**Output structure:** Same as Phase 1 (`gs.out`, `gss.out`, `ggss.out`, `p.out`, manifest, report) but with:
- `input_evidence_qc.tsv` instead of `input_gene_qc.tsv` for non-positive-controls modes
- Evidence-specific interpretation notes in `report.txt`

---

### Phase 4 — GWAS Summary Statistics

**Purpose:** Run PIGEAN with SNP-level GWAS summary statistics. The engine performs SNP-to-gene mapping internally using genome-build-matched reference files.

**Input required:**
- A GWAS summary statistics file (tab/whitespace-separated, with header)
- Explicit `--genome-build` declaration (mandatory for GWAS)

#### Basic GWAS (standard column names)

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_sumstats.tsv \
    --output runs/gwas_basic
```

> **`--genome-build` is required for GWAS.** Unlike other analysis types, GWAS will NOT default to hg19. This prevents silent coordinate-system mismatches.

#### GWAS with Custom Column Names

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build GRCh37 \
    --input ukbb_gwas.tsv \
    --gwas-chrom-col CHROM \
    --gwas-pos-col POS \
    --gwas-p-col PVAL \
    --gwas-beta-col EFFECT \
    --gwas-se-col STDERR \
    --output runs/gwas_custom_cols
```

Accepted genome build aliases (case-insensitive): `hg19`, `GRCh37`, `NCBI37` → "hg19"; `hg38`, `GRCh38` → "hg38".

#### GWAS with Global Sample Size

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_no_n.tsv \
    --gwas-n 50000 \
    --output runs/gwas_n
```

#### GWAS with Row Filter

```bash
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_multiancestry.tsv \
    --gwas-filter-col ANCESTRY \
    --gwas-filter-value EUR \
    --output runs/european
```

#### All GWAS Column Override Flags

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
| `--gwas-locus-col` | Compound chr:pos locus column |
| `--gwas-filter-col` | Row-filter column name |
| `--gwas-filter-value` | Row-filter match value |

**Output structure:** Same as Phase 1, with:
- `input_evidence_qc.tsv` (GWAS QC notes — gene QC not applicable for SNP-level input)
- GWAS-specific interpretation notes in `report.txt`
- `gwas_column_overrides` section in `run_manifest.json` (if overrides specified)

#### Error cases

```bash
# Missing --genome-build (GWAS requires it):
python3 run_pigean.py \
    --analysis gwas \
    --input gwas.tsv \
    --output /tmp/err_gwas1
# → ERROR: GWAS analysis requires --genome-build

# hg38 (reference files don't exist yet):
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg38 \
    --input gwas.tsv \
    --output /tmp/err_gwas2
# → ERROR: Required reference file not found

# Wrong column name:
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas.tsv \
    --gwas-chrom-col NONEXISTENT \
    --output /tmp/err_gwas3
# → ERROR: Specified column 'NONEXISTENT' not found in header
```

---

### Run Unit Tests

```bash
# All tests (Phase 1–4):
python3 -m pytest tests/ -v

# Individual test files:
python3 -m pytest tests/test_config.py -v
python3 -m pytest tests/test_validation.py -v
python3 -m pytest tests/test_engine_command.py -v
python3 -m pytest tests/test_gene_sets.py -v
python3 -m pytest tests/test_evidence_adapters.py -v
python3 -m pytest tests/test_gwas_adapter.py -v

# Or with unittest:
python3 -m unittest discover tests/ -v
```

**Expected result:** `173 passed`

---

### Error Cases (for verification)

These should all fail with clear error messages:

```bash
# Custom profile without files:
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets custom \
    --output /tmp/err1

# Files without format:
python3 run_pigean.py \
    --input ex/gene_list \
    --custom-gene-set-files /tmp/my_pathways.txt \
    --custom-gene-set-action replace \
    --output /tmp/err2

# Files without action:
python3 run_pigean.py \
    --input ex/gene_list \
    --custom-gene-set-files /tmp/my_pathways.txt \
    --custom-gene-set-format engine \
    --output /tmp/err3

# Nonexistent background:
python3 run_pigean.py \
    --input ex/gene_list \
    --background /does/not/exist.txt \
    --output /tmp/err4

# Nonexistent input:
python3 run_pigean.py \
    --input /does/not/exist.txt \
    --output /tmp/err5

# All fake genes (should FAIL validation):
echo -e "FAKEGENE1\nFAKEGENE2" > /tmp/fake_genes.txt
python3 run_pigean.py \
    --input /tmp/fake_genes.txt \
    --output /tmp/err6
```

---

### Quick Smoke Test (all phases, ~30 seconds)

```bash
# 1. Unit tests
python3 -m pytest tests/ -v

# 2. Golden regression test
bash tests/golden/run_golden_test.sh

# 3. Wrapper default run + golden validation
python3 run_pigean.py --input ex/gene_list --output /tmp/smoke_test --overwrite
bash tests/golden/run_golden_test.sh /tmp/smoke_test
```

---

## Changelog — Implementation Record

### Phase 0 — Golden Test and Specification (Completed 2026-09-14)

**Goal:** Freeze the current known-good run as a regression baseline.

**What was done:**
1. Ran `priors.py` 6 times (run_01 through run_06) with the same 10-gene input (`ex/gene_list`) to measure stochastic variability of the Gibbs sampler.
2. Analyzed variability across all runs:
   - Gene priors: 5–44% CV across runs (LEP most stable at 5%, RMI2 most variable at 44%)
   - `log_bf` column: perfectly deterministic across all runs (Spearman ρ = 1.0)
   - Gene set betas: very stable (pairwise ρ > 0.995)
   - Top-10 gene lists: 18–67% pairwise overlap (rank order varies due to stochasticity)
   - Overall rank stability: Spearman ρ > 0.89 for priors, > 0.95 for log_bf
3. Froze tolerance-based checkpoints in `tests/golden/checkpoints.json` (not byte-for-byte):
   - Structural: gs_data_rows=42216 (exact), gss_data_rows=41627 (exact), p_data_rows=35–36
   - 10 input gene priors with 3× observed range tolerances
   - Rank stability floors: prior min_rho=0.84736, combined min_rho=0.848117, log_bf min_rho=0.95, gene_set_beta min_rho=0.975622
   - 29 deterministic model parameters (exact match), ~6 stochastic parameters with calibrated tolerances
   - Top-20 gene/pathway overlap thresholds
4. Created `tests/golden/run_golden_test.sh` — supports both fresh engine run and validating existing outputs.
5. Created `tests/golden/expected/` with golden outputs from run_03 (representative).
6. Validated run_06 (out-of-sample, frozen after tolerances were locked): **56 PASS, 0 FAIL**.
7. Validated all 6 runs retrospectively: all pass.

**Tolerance calibration issues discovered and resolved:**
- Spearman ρ threshold 0.9 > observed min 0.897 → lowered floor to 0.84736 (with margin)
- `sigma2_cond` tolerance rounded to 0.0 by analysis script → set minimum 1e-15, later widened to 1e-9 (float precision gap of 4.1e-13)
- Top-10 overlap threshold N//2=5 > observed min 3 → changed floor to 1, expanded to top-20 with floor 8
- `p.out` row count expected 35 but 3/5 runs had 36 → added `p_data_rows_max=36` range check

**Key files created:**
- `tests/golden/checkpoints.json` — LOCKED, not to be modified without re-baselining
- `tests/golden/run_golden_test.sh`
- `tests/golden/input/positive_controls.txt`
- `tests/golden/expected/{gs,gss,ggss,p}.golden.out`
- `outputs/phase0_golden/run_01` through `run_06` — baseline runs with analysis artifacts

**Engine SHA-256 (frozen):** `da66b0c75f4554d7256daa5a2078e05113af7d746e5102d3b232ccbece2fb165`

---

### Phase 1 — Wrapper Layer + Parameterized WDL (Completed 2026-09-15)

**Goal:** Build `run_pigean.py` — a thin CLI wrapper that translates simple parameters into `priors.py` engine flags, with input QC, provenance, and a basic report. Answer the central question: *"Did the wrapper change the science?"* Answer: **No.**

**What was done:**

1. **Created `pigean/` Python package** (9 modules):
   - `pigean/__init__.py` — version string (`1.0.0`)
   - `pigean/config.py` — DEFAULTS dict (8 keys matching Phase 0), `resolve_config()` with precedence (DEFAULTS → config_file JSON → CLI args), `write_resolved_config()`
   - `pigean/references.py` — `resolve_base_dir()` (Docker vs local), `resolve_reference_paths()`, `resolve_gene_set_paths()` (mouse first, msigdb second to match Phase 0), `resolve_gene_map_path()`
   - `pigean/adapters/__init__.py` — `ValidationStatus` enum (PASS/PASS_WITH_WARNINGS/FAIL), `FileValidationResult`, `NormalizationResult`, `GeneQCResult` dataclasses, `EvidenceAdapter` base class
   - `pigean/adapters/positive_controls.py` — `PositiveControlsAdapter` with `validate_file()`, `normalize()` (preserves original in `input/`, writes normalized to `normalized/`), `build_engine_args()`
   - `pigean/validation.py` — `load_gene_map()`, `load_gene_locations()`, `count_gene_set_memberships()`, `validate_input_genes()`, `write_gene_qc_tsv()`, `compute_validation_status()`
   - `pigean/engine.py` — `build_priors_command()` (uses full flag names, e.g., `--gene-gene-set-stats-out`), `save_command()`, `run_engine()`
   - `pigean/manifest.py` — `generate_manifest()` with engine SHA-256, dependency versions, input/reference checksums, timing, three-status model; `write_manifest()`
   - `pigean/report.py` — Conservative text report: QC summary, execution status, output file row counts, interpretation notes (no "posterior probability" or "causal" claims), convergence = NOT ASSESSED

2. **Created `run_pigean.py`** (~250 lines) — thin CLI entry point orchestrating all pigean/ modules. Supports `--analysis positive-controls`, `--gene-sets default`, `--genome-build hg19`, `--preset standard`, `--config` (JSON overrides), `--overwrite`.

3. **Created unit tests** (32 tests, all pass):
   - `tests/test_config.py` — 5 tests (defaults match Phase 0, config precedence)
   - `tests/test_validation.py` — 16 tests (file validation, normalization, gene QC, combined status)
   - `tests/test_engine_command.py` — 11 tests (command structure matches Phase 0)

4. **Ran Phase 1 compatibility experiment:**
   - Original Phase 0 command → `outputs/phase1_wrapper/original_command/` → **56 PASS, 0 FAIL**
   - Phase 1 wrapper → `outputs/phase1_wrapper/wrapper_command/` → **56 PASS, 0 FAIL**
   - Command comparison documented in `outputs/phase1_wrapper/command_comparison.txt`
   - Differences are cosmetic only (full vs abbreviated flag name, absolute vs relative paths, input/ vs normalized/ with identical content)

5. **Verified failed-run behavior** (5 cases):
   - Nonexistent input → exit 1, manifest written, `FAIL`
   - Empty input → exit 1, manifest written, `FAIL`
   - All fake genes (FAKEGENE1/2/3) → exit 1, manifest written, `FAIL` ("No recognized genes")
   - Mixed valid + fake (LEP, ADIPOQ, FAKEGENE1) → `PASS_WITH_WARNINGS`, engine runs
   - Whitespace-only input → exit 1, manifest written, `FAIL`

6. **Updated infrastructure:**
   - `Dockerfile` — includes `run_pigean.py` and `pigean/`, ENTRYPOINT changed to wrapper
   - `wdl/rock_pigean.wdl` — parameterized (analysis type, preset, gene sets, genome build, memory/cpu/disk), docker tag 2.0.0
   - `README.md` — usage docs, output structure, options table, golden test instructions

**Bugs discovered and fixed:**
- `compute_validation_status()` initially returned `(status, issues)` tuple but callers expected just `status` → changed to return `ValidationStatus` only
- `generate_manifest()` API mismatch: `manifest.py` expected `file_result`/`norm_result`/`base_dir` kwargs, `run_pigean.py` passed `input_path`/`normalized_path` and omitted `base_dir` → fixed `run_pigean.py` to pass correct kwargs
- Timestamps: `run_pigean.py` used `time.time()` floats, `manifest.py` called `.isoformat()` → changed to `datetime.now()`
- `ValidationStatus` enum not JSON-serializable → added `.value` conversion in manifest.py
- `report.py` showed `✗ report.txt` in FILE LOCATIONS (checked before writing itself) → special-cased to always show `✓`

**Key artifacts produced:**
- `outputs/phase1_wrapper/wrapper_command/` — complete wrapper output (12 files)
- `outputs/phase1_wrapper/original_command/` — original Phase 0 command output for comparison
- `outputs/phase1_wrapper/command_comparison.txt` — side-by-side command diff
- `outputs/phase1_wrapper/summary.txt` — Phase 1 completion summary

**Non-negotiable constraints honored:**
- `priors.py` unchanged (SHA-256 verified: `da66b0c...`)
- `tests/golden/checkpoints.json` unchanged (Phase 0 tolerances preserved)
- No Phase 2 features implemented (no custom gene sets, no background genes, no hg38, no GMT conversion)
- Only positive-controls analysis, default gene sets, hg19

---

## Phase 2 Changelog — Custom Gene Sets + Background Genes

**Date:** 2026-09-15
**Status:** COMPLETE
**Tests:** 81 unit tests PASS, 56/56 golden checks PASS

### Capabilities Added

1. **Built-in gene-set profiles** — `--gene-sets mouse-only` and `--gene-sets msigdb-only` select subsets of the existing bundled annotation files.
2. **Custom gene sets** — User-supplied gene-set files (`--custom-gene-set-files`) in engine or GMT format (`--custom-gene-set-format`), with replace or supplement semantics (`--custom-gene-set-action`).
3. **Background genes** — `--background` specifies an eligible gene population file; maps to the engine's `--positive-controls-all-in` flag.

### Files Created

| File | Lines | Purpose |
|------|-------|---------|
| `pigean/adapters/gene_sets.py` | ~230 | GMT conversion, gene-set file validation, `prepare_custom_gene_sets()` with full legal-combinations enforcement |
| `tests/test_gene_sets.py` | ~370 | 36 tests: GMT conversion, file validation, legal combinations, resolve paths with custom, background validation |

### Files Modified

| File | Summary |
|------|---------|
| `pigean/adapters/__init__.py` | Added `GeneSetPreparationResult` and `BackgroundQCResult` dataclasses |
| `pigean/adapters/positive_controls.py` | Added `input_label`/`normalized_label` params to `normalize()` for background reuse |
| `pigean/references.py` | Added `mouse-only`/`msigdb-only` profiles to `ANNOTATION_PATHS`; extended `resolve_gene_set_paths()` with `custom_normalized_paths` and `custom_action` params |
| `pigean/config.py` | Added `custom_gene_set_format`, `custom_gene_set_action` to `CONFIG_FILE_KEYS` and `CLI_KEYS` |
| `pigean/validation.py` | Added `validate_background()` function; extended `compute_validation_status()` with `background_file_result`, `background_qc_result`, `custom_gs_prep_result` params |
| `pigean/engine.py` | Added `background_normalized_path` param to `build_priors_command()` → emits `--positive-controls-all-in` flag |
| `pigean/manifest.py` | Added `custom_gene_sets` and `background` provenance sections with checksums |
| `pigean/report.py` | Added "CUSTOM GENE SETS" and "BACKGROUND GENES" report sections |
| `run_pigean.py` | Added 4 CLI args, wired background validation/normalization, custom gene-set preparation, updated all downstream calls |
| `wdl/rock_pigean.wdl` | Added `custom_gene_set_files`, `custom_gene_set_format`, `custom_gene_set_action`, `background_gene_list` inputs with conditional flag generation |
| `tests/test_config.py` | Added `TestPhase2Config` class (5 tests for new config keys) |
| `tests/test_engine_command.py` | Added `TestBackgroundCommand` class (3 tests for `--positive-controls-all-in` flag) |
| `tests/test_validation.py` | Added `TestCombinedStatusPhase2` class (5 tests for extended combined status) |

### Key Design Decisions

1. **Legal combinations enforced in one function** — `prepare_custom_gene_sets()` is the single gatekeeper for all custom gene-set validation. Invalid combinations (e.g., files without format, format without files, custom profile without files) fail with clear error messages.
2. **Background reuses PositiveControlsAdapter** — Background gene files follow the same validate → normalize pipeline as input genes, using configurable filename labels (`input_label="background_original"`, `normalized_label="background"`).
3. **File paths stay out of config dict** — `--background` and `--custom-gene-set-files` are CLI-only file paths, not configuration settings. Only `custom_gene_set_format` and `custom_gene_set_action` enter the config precedence chain.
4. **Traceability preserved** — Custom gene-set originals land in `input/`, converted/copied versions in `normalized/`. Background files follow the same pattern. Manifest records SHA-256 checksums for all custom files.
5. **Backward compatibility** — All new function parameters default to `None`. Existing callers (Phase 1 default workflow) are completely unaffected. Golden test confirms: 56/56 PASS.

### Integration Test Results

| Scenario | Result |
|----------|--------|
| `--gene-sets custom` without `--custom-gene-set-files` | FAIL with clear message |
| `--custom-gene-set-files` without `--custom-gene-set-format` | FAIL with clear message |
| `--custom-gene-set-files` without `--custom-gene-set-action` | FAIL with clear message |
| Nonexistent `--background` file | FAIL with clear message |
| `--gene-sets mouse-only` | 1 `--X-in` flag (mouse only) |
| Custom engine-format + `--custom-gene-set-action replace` | Only custom file in `--X-in` |
| GMT file + `--custom-gene-set-action supplement` | Description stripped, 3 `--X-in` flags |
| `--background` with real gene list | `--positive-controls-all-in` in command, QC reports PLINK1 missing from background |
| Default run (golden test) | 56/56 PASS — no regression |

### Non-negotiable constraints honored
- `priors.py` unchanged (SHA-256 verified: `da66b0c...`)
- `tests/golden/checkpoints.json` unchanged (Phase 0 tolerances preserved)
- No phase-conditional logic (`if phase == 2`) in production code
- No Phase 3+ features implemented

---

### Phase 3 — Multi-Evidence Input Adapters (Completed 2026-09-21)

**Goal:** Accept precomputed gene-level evidence — Bayes factors, Z-scores, percentiles, exome associations — through a unified adapter interface, while preserving 100% backward compatibility with the positive-controls workflow.

**What was done:**

1. **Engine contract verification** — Traced all four evidence types through `priors.py` source code:
   - `_read_gene_bfs()` (line 11606): expects `Gene`/`log_bf` columns, natural-log BF scale, last-value-wins for duplicates
   - `_read_gene_zs()` (line 11730): REQUIRES `--gene-zs-id-col` and `--gene-zs-value-col`, treats values as log-odds (NOT statistical Z-scores)
   - `_read_gene_percentiles()` (line 11843): REQUIRES id/value cols, uses inverse normal transformation
   - `calculate_huge_scores_exomes()` (line 4930): auto-detects columns, requires 2 of 3 (p, beta, SE)
   - Documented in `outputs/phase3_evidence/evidence_contract.md`

2. **Created adapter framework** (`pigean/adapters/base.py`):
   - `BaseEvidenceAdapter` abstract class with `ANALYSIS_TYPE`, `ENGINE_FLAG`, `EVIDENCE_DESCRIPTION`
   - `EvidenceNormalizationResult` dataclass (original_path, normalized_path, row counts, transformations, warnings)
   - `validate_tabular_file()` shared utility for header/column checking

3. **Created four new adapters:**
   - `pigean/adapters/bayes_factor.py` — BayesFactorAdapter: validates Gene/log_bf columns (customizable), removes NA rows, warns on duplicates
   - `pigean/adapters/zscore.py` — ZScoreAdapter: REQUIRES gene_column and score_column params, documents log-odds semantics
   - `pigean/adapters/percentile.py` — PercentileAdapter: REQUIRES gene_column and score_column, supports `higher_is_better` flag
   - `pigean/adapters/exome.py` — ExomeAdapter: minimal validation (engine auto-detects), supports optional column overrides

4. **Updated adapter registry** (`pigean/adapters/__init__.py`):
   - `get_adapter(analysis_type)` function with registry dict
   - `SUPPORTED_ANALYSIS_TYPES`: positive-controls, gene-bayes-factor, gene-z-score, gene-percentile, exome

5. **Updated run_pigean.py** — full adapter-based dispatch:
   - Phase 3 CLI args: `--gene-column`, `--score-column`, `--higher-is-better`, `--exomes-gene-col`, `--exomes-p-col`, `--exomes-beta-col`, `--exomes-se-col`, `--exomes-n-col`, `--exomes-n`
   - Branching: positive-controls → NormalizationResult + validate_input_genes; other types → EvidenceNormalizationResult + validate_evidence_genes
   - `_write_fail_manifest()` helper to reduce duplication

6. **Updated supporting modules:**
   - `pigean/config.py` — Added Phase 3 config keys
   - `pigean/engine.py` — **BREAKING CHANGE**: `build_priors_command` second param changed from `normalized_input_path` (str) to `evidence_engine_args` (list); also fixed standalone flag handling
   - `pigean/validation.py` — Added `validate_evidence_genes()` and `write_evidence_qc_tsv()`
   - `pigean/manifest.py` — Added `evidence_norm_result` parameter with full evidence provenance section
   - `pigean/report.py` — Added "EVIDENCE INPUT" section + analysis-type-specific interpretation notes

7. **Test suite:**
   - `tests/test_evidence_adapters.py` — 47 new tests (registry, all adapters, backward compat)
   - `tests/test_engine_command.py` — Rewritten for evidence_args interface
   - `tests/fixtures/evidence/` — 6 fixture files
   - **Total: 128 tests PASS (81 Phase 0-2 + 47 Phase 3), 0 FAIL**

8. **Integration testing:**
   - Golden regression: **56 PASS, 0 FAIL** — positive-controls backward compat confirmed
   - gene-bayes-factor (500 genes): **COMPLETED** — 42,217 gene rows in gs.out
   - exome (500 genes): **COMPLETED** — 42,217 gene rows in gs.out
   - gene-z-score (500 genes): FAILED — engine numerical error in Gibbs sampler (scipy domain error, not wrapper bug)
   - gene-percentile (500 genes): FAILED — engine bug at line 16343 (`options.top_posterior` should be `options.gene_percentiles_top_posterior`)
   - Invalid-input tests: **7/7 PASS** (wrong columns, missing required args, nonexistent file, empty file, bad column override)

9. **Updated infrastructure:**
   - `wdl/rock_pigean.wdl` — Added Phase 3 inputs (gene_column, score_column, higher_is_better, exome overrides), docker tag 3.0.0
   - `Dockerfile` — Comment updated; `ADD pigean/ pigean/` already includes new adapter files
   - `docs/evidence_inputs.md` — User-facing documentation for all 5 evidence types

### Engine Bugs Discovered (frozen priors.py — NOT wrapper bugs)

1. **Percentile mode**: `priors.py` line 16343 references `options.top_posterior` but the attribute is `options.gene_percentiles_top_posterior`. The wrapper correctly constructs the command; the engine crashes on attribute lookup.
2. **Z-score mode**: Numerical instability in the Gibbs sampler produces a scipy domain error (`scale parameter must be positive`) with certain input distributions. The wrapper correctly validates and normalizes the input.

Both issues exist in the frozen engine (SHA-256: `da66b0c...`) and cannot be fixed without modifying priors.py.

### Non-negotiable constraints honored
- `priors.py` unchanged (SHA-256 verified)
- `tests/golden/checkpoints.json` unchanged
- Golden test: 56 PASS, 0 FAIL
- Backward compatible: default positive-controls workflow unchanged
- Every engine option traced through priors.py source
- No Phase 4/5 features implemented
- No unsupported evidence modes implemented

---

### Phase 4 — GWAS + Genome-Build Reference Packs (Completed 2026-09-21)

**Goal:** Accept GWAS summary statistics as input evidence, with mandatory explicit genome-build declaration, genome-build-matched reference file resolution, and structural prevention of build mismatches. Answer the central question: *"Can a researcher provide GWAS data with an explicit genome build and have the wrapper validate, resolve, and route the correct genomic references into PIGEAN without silently mixing coordinate systems?"* Answer: **Yes.**

**What was done:**

1. **Engine contract verification** — Traced all GWAS-related functionality through `priors.py`:
   - `calculate_huge_scores_gwas()` (line 3232): 50+ parameters, reads GWAS file, performs SNP-to-gene mapping via distance-based coordinate lookup
   - `_determine_columns()` (line 11992): auto-detects GWAS column names from header
   - `_get_col()` (line 16182): resolves column name or 1-based index
   - `_clean_chrom()` (line 12600): strips "chr" prefix, returns string
   - `_read_loc_file()` (line 12531): parses 6-column gene location files (no header)
   - Evidence priority in `read_Y()` (line 1317): GWAS gets highest priority
   - Main dispatch (line 16341): passes `gene_loc_file=options.gene_loc_file_huge if options.gene_loc_file_huge is not None else options.gene_loc_file`
   - Cataloged all 59 GWAS CLI options in the engine
   - Documented function default discrepancies: `gwas_low_p_posterior` (CLI 0.75 vs func 0.98), `gwas_high_p_posterior` (CLI 0.01 vs func 0.001), `max_clump_ld` (CLI 0.5 vs func 0.2)
   - Full contract in `outputs/phase4_gwas/gwas_engine_contract.md`

2. **Created GwasAdapter** (`pigean/adapters/gwas.py`):
   - `validate_file()`: validates tabular file, checks user-specified column overrides exist in header, supports 1-based integer column indices
   - `normalize()`: copies original to `input/gwas_original.tsv`, strips blank lines, writes to `normalized/gwas.tsv`. No data transformation — engine handles all computation
   - `build_engine_args()`: returns `["--gwas-in", path]` + all column override flags. Does NOT include gene-loc-file flags (those come from engine.py via reference_paths)
   - `COLUMN_FLAGS` dict: maps 10 wrapper param names to engine CLI flags
   - `SCALAR_FLAGS` dict: maps `gwas_n` to `--gwas-n`

3. **Genome build normalization** (`pigean/references.py`):
   - `GENOME_BUILD_ALIASES` dict: hg19/grch37/ncbi37 → "hg19", hg38/grch38 → "hg38" (case-insensitive)
   - `normalize_genome_build()`: normalizes user input to canonical form, raises ValueError for unknown builds
   - `REFERENCE_PATHS` extended with hg38 placeholder entries (files don't exist yet)
   - hg38 request produces clear `FileNotFoundError` with the missing file path

4. **Explicit genome build for GWAS** (`run_pigean.py`):
   - GWAS mode requires `--genome-build` — refusing to default to hg19 for coordinate-sensitive analysis
   - Non-GWAS modes continue to default to hg19 for backward compatibility
   - Build normalization in Step 2a (before reference resolution)
   - GWAS QC branch in Step 8: creates `EvidenceQCResult` with PASS_WITH_WARNINGS explaining that gene QC is not applicable for SNP-level input

5. **Updated adapter registry** (`pigean/adapters/__init__.py`):
   - Added `"gwas": GwasAdapter` to registry dict
   - Added `"gwas"` to `SUPPORTED_ANALYSIS_TYPES`

6. **Updated supporting modules:**
   - `pigean/config.py` — Added 11 GWAS config keys to `CONFIG_FILE_KEYS` and `CLI_KEYS`
   - `pigean/engine.py` — No changes needed: existing architecture handles GWAS args through adapter pattern; gene-loc-file flags already passed from reference_paths
   - `pigean/report.py` — Added GWAS report section (column overrides, build info, interpretation notes); also fixed QC display for all evidence types (report previously assumed `GeneQCResult` attributes on `EvidenceQCResult` objects — now handles both types correctly)
   - `pigean/manifest.py` — Handle both `GeneQCResult` and `EvidenceQCResult` in gene_qc section; added `gwas_column_overrides` provenance section

7. **Test suite:**
   - `tests/test_gwas_adapter.py` — 45 new tests across 9 test classes:
     - `TestGwasAdapterRegistry` (4 tests): registry, analysis type, engine flag
     - `TestGwasValidateFile` (11 tests): valid/invalid files, column overrides, integer indices, empty/header-only files
     - `TestGwasNormalize` (6 tests): preserves original, row counts, blank removal, content unchanged
     - `TestGwasBuildEngineArgs` (6 tests): minimal args, all overrides, scalar flags, no gene-loc flags
     - `TestGenomeBuildNormalization` (10 tests): all aliases, case insensitive, whitespace, unknown builds
     - `TestHg38ReferencesNotYetAvailable` (2 tests): hg38 fails cleanly, hg19 still works
     - `TestGwasEngineCommand` (3 tests): GWAS in full command, gene-loc files present, column overrides
     - `TestGwasConfigKeys` (3 tests): config keys present, config resolution
   - 5 test fixtures in `tests/fixtures/gwas/`: valid_gwas_hg19.tsv, gwas_custom_cols.tsv, gwas_missing_cols.tsv, gwas_no_header.txt, gwas_with_blanks.tsv
   - **Total: 173 tests PASS (128 Phase 0-3 + 45 Phase 4), 0 FAIL**

8. **Golden regression:** **56 PASS, 0 FAIL** — positive-controls backward compat confirmed.

9. **Updated infrastructure:**
   - `wdl/rock_pigean.wdl` — Added 11 GWAS inputs (`gwas_chrom_col` through `gwas_filter_value`), updated docker tag to 4.0.0, made `gene_qc` and `evidence_qc` outputs optional (`File?`)
   - `Dockerfile` — Updated comment to Phase 1-4
   - `docs/gwas_input.md` — User documentation: format requirements, genome build, column detection, examples, WDL usage, troubleshooting

### Key Design Decisions

1. **Explicit genome build, never inferred** — GWAS mode requires `--genome-build`. The engine does NOT validate coordinate system consistency. A build mismatch (hg19 GWAS + hg38 references) silently produces scientifically wrong SNP-to-gene mapping. The wrapper prevents this structurally.
2. **hg38 structural readiness without fake data** — `REFERENCE_PATHS` includes hg38 entries for when files become available, but no placeholder files were created. hg38 requests fail with `FileNotFoundError`.
3. **Minimal normalization** — The wrapper copies the GWAS file and strips blank lines but does NOT transform data values. The engine handles p/beta/SE → Z → BF → HuGE score computation and SNP-to-gene mapping internally.
4. **Column auto-detection deferred to engine** — The wrapper validates only user-specified column overrides. When no overrides are given, the engine's `_determine_columns()` handles detection. Avoids duplicating engine logic.
5. **Gene QC not applicable** — GWAS is SNP-level. The report and manifest correctly indicate gene QC doesn't apply. An `EvidenceQCResult` with `PASS_WITH_WARNINGS` is created with an explanatory note.
6. **Report/manifest QC type handling fixed** — Discovered and fixed a latent bug where `report.py` and `manifest.py` would crash for all Phase 3 evidence modes because they assumed `GeneQCResult` attributes on `EvidenceQCResult` objects. Both now handle both types correctly via attribute detection.

### Existing hg19 Reference Data (Verified)

| File | Lines | Chromosome format |
|------|-------|-------------------|
| `data/NCBI37.3.plink.gene.loc` | 19,427 | Numeric: 1-24 (23=X, 24=Y) |
| `data/refGene_hg19_TSS.subset.loc` | 34,172 | Numeric: 1-22 + X, Y |
| `data/NCBI37.3.plink.gene.exons.loc` | 336,788 | Numeric + X, Y, MT |
| `data/portal_gencode.gene.map` | 221,455 | Build-agnostic (gene symbol ↔ ID) |

### Non-negotiable constraints honored
- `priors.py` unchanged (SHA-256 verified: `da66b0c...`)
- `tests/golden/checkpoints.json` unchanged
- Golden test: 56 PASS, 0 FAIL
- Backward compatible: default positive-controls workflow unchanged
- All GWAS engine options traced through priors.py source
- Genome build never inferred from coordinates
- No silent coordinate mixing
- No liftover implemented
- hg38 NOT assumed to exist
- Phase 5 NOT implemented

---

## Quick-Start Commands

All commands run from the `rock-pigean` directory:

```bash
cd /fsx/home/l112724/gwas_pipeline/rock-pigean
```

### Phase 0 — Golden Regression Test

```bash
bash tests/golden/run_golden_test.sh
```

Runs `priors.py` with the frozen 10-gene input and validates against `tests/golden/checkpoints.json`. Expected: **56 PASS, 0 FAIL**. Add `--keep` to keep the output directory for inspection.

### Phase 1 — Wrapper (Default Settings)

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --output runs/phase1_run
```

Runs the full pipeline through `run_pigean.py` with defaults (positive-controls analysis, hg19, default gene sets, standard preset). Produces `gs.out`, `gss.out`, `ggss.out`, `p.out`, QC table, manifest, and report under `runs/phase1_run/`.

### Phase 2 — Custom Gene Sets

```bash
python3 run_pigean.py \
    --input ex/gene_list \
    --gene-sets mouse-only \
    --output runs/phase2_mouse_only
```

Runs with only the mouse phenotype gene sets (1 database instead of 2). Other Phase 2 variants:

| Variant | Flag changes |
|---|---|
| MSigDB only | `--gene-sets msigdb-only` |
| Custom file (replace defaults) | `--gene-sets custom --custom-gene-set-files /path/to/file.txt --custom-gene-set-format engine --custom-gene-set-action replace` |
| Custom file (supplement defaults) | `--gene-sets default --custom-gene-set-files /path/to/file.txt --custom-gene-set-format engine --custom-gene-set-action supplement` |
| GMT format | Use `--custom-gene-set-format gmt` instead of `engine` |

### Phase 3 — Multi-Evidence Input Types

```bash
# Gene-level Bayes factors
python3 run_pigean.py \
    --analysis gene-bayes-factor \
    --input evidence.tsv \
    --output runs/bf_run

# Gene-level scores (Z-score/log-odds mode)
python3 run_pigean.py \
    --analysis gene-z-score \
    --input scores.tsv \
    --gene-column Gene \
    --score-column Score \
    --output runs/zscore_run

# Gene-level percentiles
python3 run_pigean.py \
    --analysis gene-percentile \
    --input ranks.tsv \
    --gene-column Gene \
    --score-column rank_score \
    --output runs/percentile_run

# Gene-level percentiles (higher = better)
python3 run_pigean.py \
    --analysis gene-percentile \
    --input enrichment.tsv \
    --gene-column Gene \
    --score-column enrichment_score \
    --higher-is-better \
    --output runs/percentile_higher

# Exome associations
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --output runs/exome_run

# Exome with column overrides
python3 run_pigean.py \
    --analysis exome \
    --input exome_results.tsv \
    --exomes-gene-col GENE_SYMBOL \
    --exomes-p-col PVALUE \
    --exomes-n 50000 \
    --output runs/exome_custom
```

### Phase 4 — GWAS Summary Statistics

```bash
# Basic GWAS (engine auto-detects standard column names)
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_sumstats.tsv \
    --output runs/gwas_basic

# GWAS with custom column names
python3 run_pigean.py \
    --analysis gwas \
    --genome-build GRCh37 \
    --input ukbb_gwas.tsv \
    --gwas-chrom-col CHROM \
    --gwas-pos-col POS \
    --gwas-p-col PVAL \
    --output runs/gwas_custom

# GWAS with global sample size
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas.tsv \
    --gwas-n 50000 \
    --output runs/gwas_with_n

# GWAS with row filter
python3 run_pigean.py \
    --analysis gwas \
    --genome-build hg19 \
    --input gwas_multi.tsv \
    --gwas-filter-col ANCESTRY \
    --gwas-filter-value EUR \
    --output runs/gwas_eur
```

> **Note:** `--genome-build` is mandatory for GWAS mode. hg38 is not yet supported (reference files not available).




Phase 3 — Multi-Evidence Input Support
Add more input types — PIGEAN should accept more than just a gene list/positive controls.
Bayes Factor input — allow researchers to provide gene-level BF evidence.
Z-score input — support gene-level Z-scores.
Exome evidence — support gene-level rare-variant/exome association results.
Percentile/ranking evidence — allow ranked or percentile-based gene evidence where supported.
Create an adapter for each evidence type — each adapter validates and translates researcher input into what priors.py expects.
Standardize everything before PIGEAN — different input types become a consistent internal representation.
Validate scales carefully — for example BF may be raw BF, ln(BF), or log10(BF); never guess.
Gene-symbol/alias QC — identify unresolved or ambiguous gene names rather than silently changing them.
Regression testing — positive-control mode must still behave exactly as Phase 1/2, while each new evidence adapter gets its own tests.
Simple summary:
Phase 3 = Let researchers bring different types of gene-level evidence.


Phase 4 — GWAS + Genome Build Reference Packs
Add GWAS input support — move from gene-level inputs to SNP/variant-level summary statistics.
Create a GWAS adapter — validate and normalize GWAS summary-statistics files.
Require genome build — researcher must explicitly specify hg19/GRCh37 or hg38/GRCh38.
Validate GWAS columns — chromosome, position, P-value and other required fields such as BETA/SE/N when applicable.
Check chromosome/position formatting — catch malformed coordinates before running PIGEAN.
Create versioned reference packs — organize gene locations, TSS, exons, mappings, etc. by genome build.
Add hg38 support — current workflow is mainly hg19; Phase 4 makes both builds explicit.
Prevent build mismatch — hg38 GWAS + hg19 references should fail rather than silently produce wrong mappings.
Validate SNP-to-gene mapping behavior — confirm how variants are assigned to genes and how LD/credible-set information is handled where relevant.
Record complete reference provenance — build, reference version, checksums and files used go into the manifest.
Simple summary:
Phase 4 = Let researchers bring GWAS data safely, with the correct genome references.


Phase 5 — Interpretation Report + Convergence Diagnostics
Turn raw outputs into a scientist-readable report — gs, gss, ggss, and p become understandable summaries.
Define every reported metric — source column, scale, transformation and meaning must be verified.
Separate input genes from new candidates — don't make known positive controls look like newly discovered genes.
Show prioritized candidate genes — rank additional genes using only verified PIGEAN metrics.
Show important pathways/gene sets — summarize the biological pathways driving the model.
Explain gene-to-pathway relationships — help scientists understand why a candidate gene is being prioritized.
Add HTML reporting — move beyond raw TSV/text files to an easier report.
Capture Gibbs traces — collect the model outputs needed to evaluate sampling behavior.
Add convergence diagnostics — determine whether the statistical model actually stabilized, rather than equating successful execution with convergence.
Clearly report warnings and limitations — for example COMPLETED but CONVERGENCE NOT ESTABLISHED, instead of overstating confidence.