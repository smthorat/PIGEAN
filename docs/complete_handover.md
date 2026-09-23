# PIGEAN Complete Knowledge Transfer and Technical Handover

**Project:** rock-pigean / PIGEAN wrapper

**Current wrapper version:** 6.0.0

**Handover scope:** Phases 0–6, including the Phase 5.5 repository reorganization

**Last updated:** 2026-09-23

This is the primary knowledge-transfer document for the repository. It is
written for a researcher, engineer, or maintainer who has not participated in
the project and needs to understand why the pipeline exists, what each phase
added, how the current code works, what scientific claims are safe, how to run
and validate it, and where future changes belong.

For copy-and-paste commands, use [How to Run PIGEAN](how_to_run.md). This
document explains the design and the reasoning behind those commands.

## Table of Contents

1. [Executive summary](#1-executive-summary)
2. [The essential mental model](#2-the-essential-mental-model)
3. [Scientific concepts and interpretation boundary](#3-scientific-concepts-and-interpretation-boundary)
4. [The frozen-engine rule](#4-the-frozen-engine-rule)
5. [Repository and component map](#5-repository-and-component-map)
6. [Current end-to-end execution flow](#6-current-end-to-end-execution-flow)
7. [Phase-by-phase handover](#7-phase-by-phase-handover)
8. [Current input and compatibility contracts](#8-current-input-and-compatibility-contracts)
9. [Configuration, references, gene sets, and backgrounds](#9-configuration-references-gene-sets-and-backgrounds)
10. [Output directory and result contracts](#10-output-directory-and-result-contracts)
11. [Stability and convergence](#11-stability-and-convergence)
12. [Local, Docker, and WDL execution](#12-local-docker-and-wdl-execution)
13. [Testing and regression protection](#13-testing-and-regression-protection)
14. [Operational runbook](#14-operational-runbook)
15. [How to extend the pipeline safely](#15-how-to-extend-the-pipeline-safely)
16. [Known limitations and technical debt](#16-known-limitations-and-technical-debt)
17. [Documentation authority and historical records](#17-documentation-authority-and-historical-records)
18. [Glossary](#18-glossary)

## 1. Executive Summary

PIGEAN stands for **Probabilistic Inference of Gene Enrichment**. The scientific
engine combines gene-level evidence with gene-set annotations to estimate:

- which genes are prioritized by the model;
- which gene sets have model-estimated effects; and
- which gene-set memberships contribute to a gene's prioritization.

The original implementation is the single file `engine/priors.py`. It is large,
scientifically sensitive, stochastic, and exposes many low-level options. The
purpose of this project is not to rewrite that engine. The project places a
validated wrapper around it so a researcher can provide a supported input and
receive an auditable run directory with:

- input validation and normalization;
- gene and evidence quality control;
- controlled selection of references and gene sets;
- an exact record of the engine command;
- raw engine outputs;
- stability assessment using the engine's own criterion;
- cautious researcher-facing reports; and
- checksummed provenance.

The normal entry point is:

```bash
python3 run_pigean.py \
  --analysis positive-controls \
  --input examples/gene_list \
  --output results/runs/example_run
```

That one command executes the wrapper from input validation through reporting.
It does **not** infer what kind of evidence a file contains. The researcher must
select the correct `--analysis`, and coordinate-based GWAS data must declare a
genome build explicitly.

### What the pipeline is

- A safety, usability, interpretation, and provenance layer around PIGEAN.
- A controlled way to translate researcher-facing concepts into engine flags.
- A reproducible run contract across local Python, Docker, and Terra/Cromwell.
- A regression-protected scientific workflow whose original standard path has
  been compared against a frozen stochastic baseline.

### What the pipeline is not

- It is not a replacement or reimplementation of `priors.py`.
- It is not deterministic; the engine uses unseeded NumPy random sampling.
- It does not establish causal genes or causal pathways.
- It does not establish formal MCMC convergence through R-hat or ESS.
- It does not currently provide hg38 reference data.
- It does not safely expose every subcommand or parameter present in the
  engine.

## 2. The Essential Mental Model

Four independent choices determine a run:

1. **Evidence type** — what data the researcher supplies (`--analysis`).
2. **Model mode** — which verified engine pathway processes it (`--mode`).
3. **Gene-set annotations** — which pathway/phenotype membership matrices are
   available to the model (`--gene-sets` plus optional custom files).
4. **Genome reference** — which coordinate system and location files are used
   (`--genome-build`).

The pipeline keeps these choices separate intentionally:

```text
Researcher evidence ──> Evidence adapter ──┐
                                           │
Model choice ─────────> Mode contract ─────┼──> Frozen PIGEAN engine
                                           │
Gene sets ────────────> Annotation paths ──┤
                                           │
Genome build ─────────> Reference paths ───┘
                                                    │
                                                    v
                                      Raw outputs + reports + provenance
```

This separation prevents two common mistakes:

- treating a data format as if it automatically selected a statistical mode;
- exposing an engine feature merely because an argument exists, without first
  validating its inputs, outputs, and interpretation.

### Evidence adapters answer “what data is this?”

Adapters validate file structure, preserve the original, normalize it for the
engine, and build evidence-specific engine arguments. Examples are a simple
positive-control gene list, gene-level natural-log Bayes factors, exome
statistics, and variant-level GWAS statistics.

### Mode contracts answer “which estimator is being run?”

The `standard` mode invokes the validated outer-Gibbs workflow. The
`naive-priors` mode uses a different engine subcommand, bypasses the outer
gene-prior Gibbs update, and therefore has different output and stability
contracts. A mode/evidence combination must be explicitly marked supported or
the wrapper stops before launching the engine.

### The wrapper owns safety; the engine owns science

The wrapper decides whether a request is supported, which files to use, and how
to describe results. The frozen engine performs the underlying statistical
calculation.

## 3. Scientific Concepts and Interpretation Boundary

This section gives enough context to interpret the pipeline without claiming
more than the verified output contract supports. The complete column-by-column
contract is in [PIGEAN Engine Output Contract](output_contract.md).

### 3.1 Observed evidence: `log_bf`

`log_bf` is the gene-level observed-evidence term on a natural-log Bayes-factor
or log-odds-like scale. Its origin depends on the analysis:

- positive controls: the engine constructs the signal from the supplied
  positive-control annotation;
- gene Bayes factors: the researcher supplies natural-log Bayes factors;
- exome: the engine derives gene evidence from association statistics;
- GWAS: the engine maps variants to genes and derives a gene-level score.

Do not assume that a `log_bf` from one analysis type is operationally
interchangeable with a `log_bf` produced by another without a separate
scientific validation.

### 3.2 Gene-set-derived prior: `prior`

The model estimates gene-set effects and projects those effects through gene-set
membership to obtain a gene-level prior log-odds contribution. In standard
mode, the reported prior comes from the validated outer-Gibbs workflow. In
naive-priors mode it is calculated directly as:

```text
prior = X_orig · (beta / scale_factor)
prior = prior - mean(prior across present and missing genes)
```

These are different estimators and must not be treated as equivalent simply
because both output columns are named `prior`.

### 3.3 Combined gene evidence

The verified relationship is:

```text
combined = prior + log_bf
```

In standard mode, `combined_D` is:

```text
combined_D = exp(combined) / (1 + exp(combined))
```

The report describes `combined_D` as a **model-derived posterior probability**.
It is not an externally calibrated “true probability,” and it is not proof of
causality. Naive-priors mode does not emit `combined_D`; its report ranks genes
using `combined`.

### 3.4 Gene-set results

The main verified gene-set metrics are:

- `beta`: posterior mean gene-set effect corrected for correlation between
  gene sets; the primary multivariate effect estimate;
- `avg_postp`: posterior inclusion probability for a non-zero gene-set effect;
- `P`: marginal regression p-value from the initial filter, not from the
  multivariate Gibbs estimate; and
- `N`: gene-set size.

`P` and `avg_postp` answer different questions and should not be expected to
rank gene sets identically.

### 3.5 Safe language

Reports intentionally use words such as **prioritized**, **model-derived**,
**candidate**, **associated**, and **enriched**. They avoid claims such as
**causal**, **confirmed**, **validated gene**, or **true probability**.

For positive-control runs, input genes are part of the training signal. Their
high scores are expected and are not independent discoveries. The report splits
them from non-input candidate genes for this reason.

## 4. The Frozen-Engine Rule

`engine/priors.py` is the scientific engine and is treated as frozen.

```text
Frozen engine SHA-256:
da66b0c75f4554d7256daa5a2078e05113af7d746e5102d3b232ccbece2fb165
```

### Why it is frozen

The wrapper was developed to generalize and make the existing scientific
workflow auditable without silently changing its results. Modifying the engine
would mix two kinds of work:

- engineering changes to validation, orchestration, and reporting; and
- scientific changes to the statistical implementation.

Keeping the boundary strict lets the golden test answer the key regression
question: “Did the wrapper change the established standard workflow?”

### What to do when an engine defect is found

1. Trace and document the defect precisely.
2. Mark the affected wrapper combination `ENGINE_BLOCKED` or `NOT_VERIFIED`.
3. Reject it before engine execution.
4. Do not patch `priors.py` as an incidental wrapper fix.
5. If a scientific engine change is approved later, create a separate engine
   version, establish a new baseline, and document the new checksum and output
   contract.

This rule is why Z-score and percentile adapters exist in the code but are not
currently runnable through either exposed model mode.

## 5. Repository and Component Map

### 5.1 Top-level layout

| Path | Responsibility | Change policy |
|---|---|---|
| `run_pigean.py` | User-facing CLI and end-to-end orchestration | Change when workflow order or exposed interface changes |
| `pigean/` | Wrapper implementation | Normal development area |
| `engine/priors.py` | Frozen scientific engine | Do not modify without a separately governed scientific rebaseline |
| `data/` | Bundled gene maps, hg19 locations, and built-in gene sets | Treat as versioned scientific inputs; checksum changes are material |
| `tests/` | Unit, fixture, and golden regression coverage | Update with every contract change |
| `workflows/` | Terra/Cromwell WDL interface | Keep synchronized with CLI flags and container version |
| `Dockerfile` | Reproducible runtime image | Keep engine, wrapper, data, and version tag aligned |
| `docs/` | Current contracts, guides, and historical phase records | Update whenever behavior or interpretation changes |
| `results/` | Generated runs and development artifacts | Ignored by Git; never use as source code |

### 5.2 Wrapper modules

| Module | Why it exists | What it does |
|---|---|---|
| `pigean/config.py` | One reproducible settings resolution path | Merges defaults, preset, JSON config, and explicit CLI values; writes `resolved_config.json` |
| `pigean/references.py` | Prevent scattered or environment-specific path logic | Locates repo vs `/app`, normalizes genome-build aliases, verifies reference and gene-set files |
| `pigean/adapters/base.py` | Standardize evidence integrations | Defines tabular validation, normalization/QC result types, and the adapter interface |
| `pigean/adapters/positive_controls.py` | Preserve the original gene-list workflow | Validates one-gene-per-line text, removes blanks/duplicates, builds `--positive-controls-in` |
| `pigean/adapters/bayes_factor.py` | Safely accept gene-level evidence | Requires gene and natural-log-BF columns, filters invalid values, normalizes headers |
| `pigean/adapters/exome.py` | Route rare-variant association evidence | Preserves tabular input and passes optional gene/p/beta/SE/N mappings |
| `pigean/adapters/gwas.py` | Route SNP-level summary statistics | Validates headers/overrides, preserves values, builds GWAS-specific arguments |
| `pigean/adapters/zscore.py` | Preserve a traced but blocked engine contract | Validates/normalizes input; mode layer prevents production execution |
| `pigean/adapters/percentile.py` | Preserve a traced but blocked engine contract | Validates/normalizes input; mode layer prevents production execution |
| `pigean/adapters/gene_sets.py` | Make custom annotations auditable | Validates custom files, converts GMT, preserves originals, enforces legal combinations |
| `pigean/validation.py` | Stop bad inputs before expensive computation | Performs gene-map, coordinate, annotation, background, and combined-status QC |
| `pigean/modes/` | Separate statistical pathway from evidence type | Registers verified modes, compatibility, stability applicability, and subcommands |
| `pigean/engine.py` | Centralize the exact engine invocation | Builds the command list, saves `priors_command.txt`, runs the subprocess, captures logs |
| `pigean/parsers.py` | Turn engine TSVs into typed structures | Parses and indexes `gs.out`, `gss.out`, `ggss.out`, and `p.out` without pandas |
| `pigean/convergence.py` | Use precise stability language | Parses parameters/logs and writes the engine-specific stability assessment |
| `pigean/interpretation.py` | Enforce a verified scientist-facing contract | Ranks genes/gene sets, links genes to pathways, emits caveats, avoids forbidden claims |
| `pigean/report.py` | Make results consumable without hiding raw data | Generates `report.txt` and self-contained `report.html` |
| `pigean/manifest.py` | Make every completed/handled run auditable | Records versions, hashes, settings, validation, execution, stability, and summaries |

### 5.3 Reference data

| File | Purpose |
|---|---|
| `data/portal_gencode.gene.map` | Maps source identifiers to gene symbols for recognition/QC |
| `data/NCBI37.3.plink.gene.loc` | hg19/GRCh37 gene coordinates |
| `data/refGene_hg19_TSS.subset.loc` | hg19 transcription-start-site coordinates used in HuGE mapping |
| `data/NCBI37.3.plink.gene.exons.loc` | hg19 exon coordinates |
| `data/gene_set_list_mouse_2024.txt` | Mouse phenotype annotation gene sets |
| `data/gene_set_list_msigdb_nohp.txt` | MSigDB annotation gene sets |

## 6. Current End-to-End Execution Flow

The orchestration lives in `run_pigean.py`. The order is deliberate because an
expensive scientific run should begin only after configuration, compatibility,
files, references, and QC are resolved.

### Step 1 — Prepare the output directory

The wrapper converts `--output` to an absolute path. A nonempty directory is
rejected unless `--overwrite` is present. With `--overwrite`, the directory is
deleted and recreated, so callers must verify the target path carefully.

### Step 2 — Resolve configuration

`pigean.config.resolve_config()` applies this precedence, from lowest to
highest:

```text
hard-coded defaults -> preset -> JSON config -> explicit CLI values
```

Only recognized keys are copied from JSON or CLI into the resolved scientific
configuration. Input/output/background/custom-file paths stay on the CLI.

### Step 3 — Normalize and validate the genome build

Aliases are normalized case-insensitively:

- `hg19`, `GRCh37`, `NCBI37` -> `hg19`
- `hg38`, `GRCh38` -> `hg38`

GWAS requires an explicit `--genome-build`; the wrapper will not silently use
the default because a coordinate mismatch can produce scientifically incorrect
SNP-to-gene mapping. hg38 aliases are recognized structurally, but the required
files are not bundled, so hg38 resolution fails clearly.

### Step 4 — Write the resolved configuration

`resolved_config.json` records the effective values after precedence and build
normalization. This is the configuration to cite when reproducing a run, not
the user's partial JSON alone.

### Step 5 — Resolve and validate the model mode

`pigean.modes.get_mode()` rejects unknown modes. The selected mode then checks
the evidence/mode compatibility matrix and rejects hidden internal parameters
such as `anchor`, `phi`, and `alpha0`. There is no silent fallback to standard
mode.

### Step 6 — Select the evidence adapter

`pigean.adapters.get_adapter()` maps `--analysis` to an adapter instance. The
adapter owns file-level validation, normalization, and evidence-specific engine
arguments.

### Step 7 — Validate and preserve the primary input

The adapter checks that the file exists, is readable, nonempty, textual, and
structurally compatible with the selected analysis. The original is copied
under `input/`; the engine-ready version is written under `normalized/`.

Normalization is intentionally conservative:

- positive controls: trim, remove blanks, and deduplicate;
- Bayes factors: remove invalid rows and normalize selected headers;
- exome/GWAS: remove blank lines but do not transform scientific values;
- custom GMT gene sets: remove the description field to create engine format.

### Step 8 — Prepare optional background and custom gene sets

A background gene list reuses the positive-control adapter and is preserved and
normalized separately. Custom gene sets are validated, copied, and optionally
converted before they can participate in reference resolution.

### Step 9 — Resolve references and annotation files

`pigean.references` selects absolute paths relative to `/app` in Docker or the
repository root locally. Every required file must exist. Gene-set ordering is
preserved because the Phase 0 command used mouse annotations before MSigDB.

### Step 10 — Run QC and combine validation status

Positive-control QC checks recognition, genomic coordinates, and annotation
membership for every input gene. Gene-level tabular evidence checks gene
recognition. GWAS is SNP-level, so the QC output explicitly says gene QC is not
applicable; SNP-to-gene mapping occurs inside the engine.

The combined status is:

- `PASS`: all checks passed;
- `PASS_WITH_WARNINGS`: execution may continue, but issues are recorded; or
- `FAIL`: the engine is not launched.

Any fatal primary-file, QC, background, or custom-gene-set result makes the
combined status fail.

### Step 11 — Build and save the exact engine command

The wrapper combines:

- Python and the frozen engine path;
- the mode's engine subcommand;
- mode-specific arguments;
- ordered `--X-in` gene-set files;
- the gene map and genome references;
- locked scientific defaults;
- output destinations;
- evidence-adapter arguments;
- optional background; and
- optional standard-mode trace output.

The exact command is written to `priors_command.txt` before execution.

### Step 12 — Execute the engine

The command is run with `subprocess.run`. Standard output and standard error are
merged into `pigean_run.log`. The wrapper records the engine exit code and does
not treat file creation alone as proof of success.

### Step 13 — Parse outputs

After a successful engine exit, parsers load the four core TSV outputs:

- genes are sorted by `combined_D`, falling back to `combined` and then `prior`;
- gene sets are sorted by absolute `beta`;
- gene-gene-set entries receive indexes in both directions; and
- repeated parameter versions are stored in order.

### Step 14 — Assess stability

Standard mode uses `p.out` plus log parsing when available. Naive-priors returns
a deliberately shaped `NOT_APPLICABLE` result because it does not run the outer
Gibbs loop being assessed.

### Step 15 — Interpret results

The interpretation layer extracts only verified fields. It separates input
positive controls from candidate genes, ranks gene sets, constructs top
gene/pathway links, selects mode-appropriate ranking, and supplies caveats.

### Step 16 — Write provenance and reports

The final `run_manifest.json` contains engine/dependency versions, checksums,
resolved analysis, mode metadata, QC, execution status, timing, stability, and
summary counts. Text and HTML reports are then generated from the structured
interpretation.

### Failure behavior

Failures after mode resolution normally write a manifest with validation
`FAIL` and execution `NOT_RUN`. Engine failures record a nonzero exit and retain
the log and any partial outputs for diagnosis. Two very early checks—an invalid
genome-build token and a GWAS request without explicit `--genome-build`—currently
exit before the normal failure-manifest helper is called. Maintainers should be
aware of this exception when relying on “every attempted run has a manifest.”

## 7. Phase-by-Phase Handover

The phases are development history, not runtime switches. Production code must
never branch on a phase number. Each phase extended the same end-to-end command
while protecting the standard baseline.

| Phase | Capability milestone | Current state |
|---|---|---|
| 0 | Stochastic golden baseline | Locked |
| 1 | Wrapper, validation, provenance, Docker/WDL | Locked |
| 2 | Custom gene sets and backgrounds | Locked |
| 3 | Multi-evidence adapters | Implemented; Z-score/percentile engine routes subsequently blocked |
| 4 | GWAS and genome-build safety | Locked for hg19/GRCh37 |
| 5 | Interpretation, HTML report, and stability terminology | Locked |
| 5.5 | Repository reorganization | Complete; no scientific change |
| 6 | Advanced-mode contracts and naive-priors | Complete with documented compatibility limits |
| 7 | Broader production infrastructure from the original roadmap | Not implemented in this repository |

### Phase 0 — Golden baseline and scientific regression specification

**Why:** The engine is stochastic, so a byte-for-byte output comparison would
fail even when behavior is correct. Before adding a wrapper, the project needed
an empirical definition of “the established workflow still behaves the same.”

**What was added:**

- six independent runs of the original positive-control command;
- a representative set of expected raw outputs;
- calibrated tolerance checkpoints;
- rank-correlation and top-N overlap checks; and
- `tests/golden/run_golden_test.sh`.

**How it works:** The golden test checks exact structural invariants, tolerant
numeric checkpoints, rank correlations, overlap of top genes/pathways, and
model parameters. It supports both a fresh engine run and validation of an
existing output directory.

**Key files:**

| File | Role |
|---|---|
| `tests/golden/checkpoints.json` | Locked structural, numeric, correlation, and overlap thresholds |
| `tests/golden/input/positive_controls.txt` | Baseline ten-gene input |
| `tests/golden/expected/` | Representative output tables |
| `tests/golden/run_golden_test.sh` | Runner and 56-check validator |

**Verification:** All six calibration runs passed the locked criteria; the
fresh regression run before the Phase 6 commit passed 56/56 checks.

**Handover rule:** Do not loosen a checkpoint merely to make a new change pass.
Changing the baseline requires a documented scientific reason and a new
multi-run calibration.

### Phase 1 — CLI wrapper, validation, provenance, and parameterized WDL

**Why:** The original engine command exposed internal details, used hardcoded
paths, and provided no unified input QC or audit trail. Researchers needed a
stable interface that reproduced the established command.

**What was added:**

- `run_pigean.py` as the public entry point;
- default/preset/config resolution;
- local/Docker reference resolution;
- positive-control validation and normalization;
- gene recognition, coordinate, and annotation QC;
- exact command construction and logging;
- a run manifest and basic report;
- Docker ENTRYPOINT integration; and
- a parameterized WDL workflow.

**How it preserves science:** The wrapper emits the same effective standard
engine settings and evidence path as the Phase 0 command. Original-command and
wrapper-command outputs both passed the golden test.

**Key components:** `config.py`, `references.py`,
`adapters/positive_controls.py`, `validation.py`, `engine.py`, `manifest.py`,
`report.py`, `run_pigean.py`, `Dockerfile`, and `workflows/rock_pigean.wdl`.

**Important design choices:**

- original and normalized inputs are both retained;
- warnings do not silently become failures, and failures do not launch the
  engine;
- the effective command is saved rather than reconstructed later; and
- scientific defaults live centrally in `DEFAULTS`.

### Phase 2 — Built-in/custom gene sets and background genes

**Why:** Researchers need to change the annotation universe or define the
eligible population without manually editing the engine command.

**What was added:**

- `mouse-only` and `msigdb-only` built-in profiles;
- custom engine-format and GMT gene-set files;
- `replace` and `supplement` behavior;
- multiple custom files;
- positive-control background genes; and
- provenance/QC for all of the above.

**How custom gene sets work:**

1. `prepare_custom_gene_sets()` enforces legal flag combinations.
2. Each source file is preserved under `input/`.
3. Engine-format files are copied; GMT files have their description column
   removed and are written under `normalized/`.
4. `resolve_gene_set_paths()` constructs the final ordered `--X-in` list.

**Replace vs supplement:**

- `replace`: use only supplied custom files;
- `supplement`: use the selected built-in profile plus supplied files;
- `--gene-sets custom --custom-gene-set-action supplement`: uses the default
  built-in profile as the base plus the custom files.

**How backgrounds work:** A background is a one-gene-per-line list passed to
the engine as `--positive-controls-all-in`. QC counts mapped/unmapped genes and
warns when positive-control input genes are missing from the background.

**Key components:** `adapters/gene_sets.py`, extensions in `references.py`,
`validation.py`, `engine.py`, `manifest.py`, `report.py`, the CLI, and WDL.

**Handover rule:** Backgrounds are a positive-control concept in the current
validated interface. Do not use `--background` with other evidence types
without first adding explicit compatibility validation and tests.

### Phase 3 — Multi-evidence adapter architecture

**Why:** The engine can ingest more than positive controls, but each format has
different columns, scale semantics, duplicate behavior, and normalization.
Embedding all format logic in the CLI would be brittle.

**What was added:**

- `BaseEvidenceAdapter` and shared tabular validation;
- Bayes-factor, Z-score, percentile, and exome adapters;
- a registry selected by `--analysis`;
- evidence-specific normalization provenance and QC; and
- corresponding CLI/WDL column mappings.

**Adapter contract:**

```text
validate_file(input, params)
normalize(input, output_dir, params)
build_engine_args(normalized_input, params)
```

**Current outcome by adapter:**

- Bayes factor: validated and supported in standard mode.
- Exome: validated and supported in standard mode.
- Z-score: wrapper adapter is implemented, but the frozen engine can produce a
  SciPy domain error from numerical instability; both modes reject it.
- Percentile: wrapper adapter is implemented, but the frozen engine references
  an incorrect option attribute; both modes reject it.

This is an important handover distinction: presence in the adapter registry or
CLI choices means “the input contract exists,” not “the complete scientific run
is supported.” The mode compatibility layer is the final authority.

### Phase 4 — GWAS evidence and genome-build protection

**Why:** GWAS is variant-level and coordinate-sensitive. Silently pairing an
hg19 GWAS file with hg38 references, or vice versa, could produce plausible but
scientifically wrong gene mapping.

**What was added:**

- `GwasAdapter` for summary-statistic inputs;
- explicit GWAS column mappings and row filtering;
- build alias normalization;
- mandatory explicit build declaration for GWAS;
- build-matched location/TSS/exon resolution; and
- GWAS-specific QC, reporting, manifest data, WDL inputs, and tests.

**How it works:** The wrapper validates the file and requested column names but
does not reproduce the engine's statistical transformations. It passes the
normalized file and resolved hg19 references to the engine, which performs
variant evidence processing and SNP-to-gene mapping.

**Why GWAS gene QC is different:** The input contains variants, not genes. The
wrapper therefore records that gene QC is not applicable; genes emerge from the
engine's coordinate-based mapping.

**Current reference boundary:** hg19/GRCh37 is operational. hg38/GRCh38 aliases
exist so the interface is structurally ready, but the files are intentionally
absent. No liftover is implemented.

**Key components:** `adapters/gwas.py`, `references.py`, CLI/config additions,
report/manifest handling, `docs/gwas_input.md`, WDL, and
`tests/test_gwas_adapter.py`.

### Phase 5 — Output parsing, interpretation, HTML reporting, and stability

**Why:** Raw engine tables are extensive and easy to misinterpret. Execution
success also does not prove that the engine's own stability criterion was met.

**What was added:**

- typed, standard-library parsers for all four raw outputs;
- a field-level output contract with VERIFIED/UNVERIFIED classifications;
- top gene and gene-set summaries;
- gene-to-pathway and pathway-to-gene links;
- cautious scientific caveats and forbidden terminology;
- a self-contained HTML report;
- PIGEAN max-fractional-SEM stability assessment; and
- a separate formal-MCMC status.

**Interpretation policy:** A value is interpreted only when its definition,
scale, direction, and anti-interpretation have been traced. Unverified values
may be left raw or omitted; they are not dressed up with a convenient label.

**Stability terminology correction:** Early development language used the word
“converged.” The audit determined that the observable calculation is the
engine's specific max-fractional-SEM criterion, not a complete R-hat/ESS
analysis. Status names were changed accordingly.

**Key components:** `parsers.py`, `convergence.py`, `interpretation.py`, major
extensions in `report.py` and `manifest.py`, `docs/output_contract.md`, and
`docs/convergence.md`.

### Phase 5.5 — Repository reorganization

**Why:** Development artifacts and production code were mixed at the root, and
the scientific engine was not visually separated from the wrapper.

**What changed:**

- `priors.py` moved to `engine/priors.py`;
- WDL moved to `workflows/`;
- example input moved to `examples/`;
- scripts moved under `scripts/`;
- generated artifacts moved under `results/`; and
- phase records moved under `docs/phases/`.

All path consumers—wrapper, Dockerfile, WDL, manifest hashing, tests, and
scripts—were updated. The engine and reference data remained byte-for-byte
unchanged.

**Handover rule:** A future reorganization must begin with a path dependency
audit. See [Path Dependency Audit](path_dependency_audit.md).

### Phase 6 — Verified advanced-mode framework and naive-priors

**Why:** The engine exposes many subcommands and internal controls, but exposing
them directly would imply they share the standard input, output, stability, and
interpretation contracts. They do not.

**What was added:**

- a high-level mode registry separate from evidence adapters;
- explicit evidence/mode compatibility states;
- the default `standard` mode;
- the verified positive-control-only `naive-priors` mode;
- mode metadata in configuration, manifest, reports, and WDL;
- mode-aware sorting and interpretation; and
- pre-engine rejection of unverified combinations and internal knobs.

**Standard mode:**

- engine subcommand: `gibbs`;
- Phase 1–5 behavior;
- supports positive controls, gene Bayes factors, exome, and GWAS;
- standard output schema and Phase 5 interpretation contract; and
- PIGEAN stability criterion applies.

**Naive-priors mode:**

- engine subcommand: `naive_priors`;
- supported only with positive controls;
- estimates gene-set effects with the inner sampler, then computes a
  mean-centered prior directly;
- bypasses the outer gene-prior Gibbs update;
- emits all four core files, but `gs.out` lacks `combined_D` and adjustment
  columns under the locked settings;
- ranks genes using `combined`; and
- reports outer-Gibbs stability and formal convergence as `NOT_APPLICABLE`.

**Investigated but not exposed:** factor/`naive_factor`, PheWAS, anchor
variants, `phi`, `alpha0`, and `beta0`. Their input matrices, masks, output
schemas, or parameter sensitivity contracts are not compatible with the
current researcher-facing model. `beta0` is parsed by the engine but is not
consumed by the active calculation, so exposing it would be especially
misleading.

**Key components:** `pigean/modes/`, mode-aware changes throughout CLI,
parsers, convergence, interpretation, reports, manifest, Docker/WDL, and
`tests/test_advanced_modes.py` plus `tests/test_naive_priors.py`.

## 8. Current Input and Compatibility Contracts

### 8.1 Authoritative compatibility matrix

| Researcher evidence | `--analysis` | `standard` | `naive-priors` |
|---|---|---:|---:|
| One gene per line | `positive-controls` | SUPPORTED | SUPPORTED |
| Gene-level natural-log Bayes factors | `gene-bayes-factor` | SUPPORTED | NOT_VERIFIED; rejected |
| Gene-level exome/rare-variant statistics | `exome` | SUPPORTED | NOT_VERIFIED; rejected |
| SNP-level GWAS summary statistics | `gwas` | SUPPORTED | NOT_VERIFIED; rejected |
| Gene-level score values | `gene-z-score` | ENGINE_BLOCKED; rejected | ENGINE_BLOCKED; rejected |
| Gene-level ranks/percentiles | `gene-percentile` | ENGINE_BLOCKED; rejected | ENGINE_BLOCKED; rejected |

The CLI lists all adapter values, including the two blocked values. Always use
this matrix or `pigean/modes/standard.py` and `naive_priors.py` to determine
whether an end-to-end run is supported.

### 8.2 Positive controls

- Format: one gene symbol per line, no header.
- Normalization: trim whitespace, remove blanks, preserve first occurrence.
- QC: gene-map recognition, coordinates, annotation membership.
- Optional background: one gene per line; should contain every input gene.
- Duplicate behavior: wrapper removes duplicates.

### 8.3 Gene Bayes factors

- Format: whitespace/tab-delimited with a header.
- Defaults: `Gene` and `log_bf`.
- Scale: natural logarithm of BF, not log10(BF).
- Custom mappings: `--gene-column`, `--score-column`.
- Invalid/NA/nonfinite rows: removed or reported during normalization.
- Duplicate behavior: engine uses the last occurrence.

### 8.4 Exome/rare-variant evidence

- Format: whitespace/tab-delimited with a header.
- Required statistical information: at least two of p-value, beta, and SE;
  sample size may support derivation.
- Standard headers are engine-detected; wrapper flags can override gene, p,
  beta, SE, and N columns or supply a global N.
- Wrapper does not transform values; the engine derives gene evidence.
- Duplicate behavior: engine keeps the first occurrence.

### 8.5 GWAS summary statistics

- Format: variant-level whitespace/tab-delimited data with a header.
- Coordinates: chromosome+position or a compound locus column.
- Statistics: enough of p-value, beta, SE, and sample size for the engine.
- Build: explicit hg19/GRCh37 required with the current data bundle.
- Column selectors: names or 1-based indices.
- Wrapper normalization: preserves values and removes blank lines.
- Gene QC: not applicable before engine mapping.

Detailed commands and all GWAS flags are in [GWAS Input](gwas_input.md).

## 9. Configuration, References, Gene Sets, and Backgrounds

### 9.1 Locked defaults

| Setting | Default | Purpose |
|---|---:|---|
| `analysis` | `positive-controls` | Evidence adapter |
| `mode` | `standard` | Statistical pathway |
| `preset` | `standard` | Named settings group; currently no alternate preset |
| `genome_build` | `hg19` | Reference pack |
| `gene_sets` | `default` | Annotation profile |
| `max_num_gene_sets` | `5000` | Engine gene-set cap |
| `gene_filter_value` | `1` | Engine gene filter |
| `gene_set_filter_value` | `0.01` | Engine gene-set filter |
| `debug_level` | `3` | Engine log detail |

These defaults reproduce the established standard command. Treat changes to
the last four scientific settings as behavior changes requiring regression and
interpretation review.

### 9.2 JSON configuration

JSON may contain recognized settings and column mappings. Explicit CLI values
override JSON. File locations are intentionally kept out of JSON:

- primary input and output;
- background file; and
- custom gene-set file paths.

This avoids ambiguity about which layer owns path localization, particularly
in WDL and Docker.

### 9.3 Gene-set profiles

| Profile | Final built-in annotation files |
|---|---|
| `default` | Mouse, then MSigDB |
| `mouse-only` | Mouse only |
| `msigdb-only` | MSigDB only |
| `custom` | Custom-only for `replace`; default+custom for `supplement` |

Custom files must always declare `--custom-gene-set-format` (`engine` or `gmt`)
and `--custom-gene-set-action` (`replace` or `supplement`).

### 9.4 Reference checksums

The run manifest hashes the engine, gene map, selected gene sets, and location
references. This matters because the same CLI settings with different data
files are not the same scientific run.

## 10. Output Directory and Result Contracts

### 10.1 Directory structure

The exact set varies by analysis and mode, but a successful run normally has:

```text
<output>/
├── input/                       preserved researcher inputs
├── normalized/                  engine-ready copies
├── resolved_config.json         effective settings
├── input_gene_qc.tsv            positive-control QC, or
├── input_evidence_qc.tsv        tabular/GWAS QC
├── priors_command.txt           exact engine command
├── pigean_run.log               engine stdout + stderr
├── gs.out                       gene statistics
├── gss.out                      gene-set statistics
├── ggss.out                     gene ↔ gene-set statistics
├── p.out                        engine parameters
├── gss_trace.out                optional standard-mode trace
├── convergence.json             stability/convergence applicability
├── run_manifest.json            provenance and statuses
├── report.txt                   text interpretation
└── report.html                  self-contained HTML interpretation
```

### 10.2 `gs.out` — gene statistics

Use this file to inspect gene-level evidence and prioritization. Key verified
columns include:

| Column | Meaning |
|---|---|
| `Gene` | Gene symbol |
| `prior` | Gene-set-derived prior log-odds contribution; estimator is mode-specific |
| `log_bf` | Observed-evidence term |
| `combined` | `prior + log_bf` |
| `combined_D` | Standard-mode sigmoid of `combined`; absent in naive-priors |
| `positive_control` | Identifies supplied training genes when present |
| `N` | Number of gene-set memberships |
| `Chrom`, `Start`, `End` | Reference coordinates |

Adjusted and internal columns are not automatically safe to interpret. Consult
the output contract before using them.

### 10.3 `gss.out` — gene-set statistics

Use this file for gene-set effects. Reports rank by absolute `beta`, not by
the marginal `P` column. `avg_postp` is a posterior inclusion probability and
is not a p-value.

### 10.4 `ggss.out` — gene/gene-set links

This sparse long-format file connects a gene with gene sets whose effects pass
the engine's output filter. It supports explanations such as “which reported
gene sets contribute to this candidate gene?” It is not a complete dense
membership matrix.

### 10.5 `p.out` — parameters and execution summary

Parameters may have multiple versions after restarts or hyperparameter
adjustment. The parser retains all versions and uses the latest where needed.
`num_gibbs_iter` is zero-indexed: 499 means 500 iterations ran.

### 10.6 Reports vs raw outputs

`report.html` is the normal researcher entry point. Raw outputs remain the
source for downstream quantitative analysis. The report intentionally selects
only verified metrics and adds caveats; it does not replace the raw tables.

### 10.7 Manifest statuses

Keep three ideas separate:

- **validation status:** was the request/input acceptable?
- **execution status:** did the engine process complete?
- **stability status:** what does the engine-specific diagnostic say?

A run can validate successfully and execute successfully while having an
uncertain or not-met stability status.

## 11. Stability and Convergence

### 11.1 What is actually assessed

The engine calculates a maximum fractional standard error of the mean across
running cross-chain gene estimates. The locked threshold is `0.01`.

The wrapper uses two evidence tiers:

- `SUMMARY_BASED`: `p.out` iteration and restart information;
- `TLOG_BASED`: log parsing provides the final max-fractional-SEM value and can
  verify the threshold directly.

### 11.2 Canonical stability statuses

- `PIGEAN_STABILITY_CRITERION_MET`
- `PIGEAN_STABILITY_CRITERION_UNCERTAIN`
- `PIGEAN_STABILITY_CRITERION_NOT_MET`
- `NOT_ASSESSED`
- `ENGINE_FAILED`
- `NOT_APPLICABLE`

Hitting the iteration cap and meeting the SEM criterion are separate facts. A
run can do both.

### 11.3 What is not assessed

Standard mode reports formal MCMC convergence as
`NOT_FORMALLY_ASSESSED`. The engine does not retain the independent chain
traces required for standard R-hat and effective-sample-size calculations.
Burn-in R-like log values may be reported as supporting information but are not
promoted to a formal convergence claim.

Naive-priors reports both outer-Gibbs stability and formal convergence as
`NOT_APPLICABLE`, while separately warning that its inner gene-set sampler is
still stochastic.

See [PIGEAN Stability and Convergence](convergence.md) for the full decision
logic.

## 12. Local, Docker, and WDL Execution

### 12.1 Local Python

Requirements are Python 3, NumPy, and SciPy. Run from the repository root so
the wrapper can resolve the engine and data relative to its package location.
The engine is compute-intensive; runtime is much longer than the unit tests.

### 12.2 Docker

The Docker image installs Python and SciPy, copies the hg19 references, frozen
engine, wrapper package, and documentation, then sets the wrapper as ENTRYPOINT.
Consequently, arguments after the image name are wrapper arguments.

Build locally with:

```bash
docker build -t rock-pigean:6.0.0 .
```

Input and output paths must be bind-mounted. The recommended full command is in
[How to Run PIGEAN](how_to_run.md).

### 12.3 WDL/Terra/Cromwell

`workflows/rock_pigean.wdl` exposes the CLI's analysis, mode, gene-set,
background, column, GWAS, and trace settings. It localizes files and runs the
same wrapper inside the configured image.

Default runtime request:

- 4 CPUs;
- 8 GB memory; and
- 50 GB local disk.

The WDL image reference is `gcr.io/nitrogenase-docker/rock-pigean:6.0.0`.
Publishing a repository commit does not automatically publish that container;
image availability and tag contents must be checked separately.

## 13. Testing and Regression Protection

### 13.1 Unit tests

Run:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
```

The Phase 6 repository has 346 tests covering configuration, validation,
command construction, gene sets, all adapters, GWAS, parsers, stability,
interpretation, modes, and naive-priors behavior.

Tests use fixtures and mocks for fast contract checks; they do not replace a
fresh engine run.

### 13.2 Golden regression

Run a fresh standard baseline:

```bash
bash tests/golden/run_golden_test.sh
```

Or validate an existing engine output directory:

```bash
bash tests/golden/run_golden_test.sh /absolute/path/to/output
```

The suite contains 56 checks. The last pre-handover fresh run passed 56/56.
Expect the fresh run to take minutes because it executes the full stochastic
engine.

### 13.3 What must be tested for every change

| Change | Minimum evidence |
|---|---|
| Documentation only | Link/path review and `git diff --check` |
| Validation/config/adapter change | Targeted tests plus full unit suite |
| Command/reference/default change | Full unit suite plus golden regression |
| Mode or interpretation change | Mode tests, schema tests, report tests, representative end-to-end run |
| Docker/WDL change | Local image build, container run, and WDL syntax/interface review |
| Engine/reference-data change | Separate scientific approval and new baseline calibration |

### 13.4 Why raw diffs are wrong for the golden test

The engine does not set `np.random.seed()`. Gene priors, gene-set effects, and
some parameters vary. The golden test therefore protects distributions,
structure, ranking, and calibrated values—not identical bytes.

## 14. Operational Runbook

### 14.1 Before starting a run

1. Confirm the scientific evidence type and its scale.
2. Confirm the selected mode is supported for that evidence.
3. Confirm genome build, especially for GWAS.
4. Inspect headers and decide whether column overrides are needed.
5. Decide built-in vs custom gene sets and replace vs supplement.
6. Use a new output directory or deliberately approve `--overwrite`.
7. For Docker, confirm all host input paths are mounted read-only and output is
   mounted read-write.

### 14.2 After a run

1. Check the process exit code and `execution.status` in the manifest.
2. Read validation warnings and the appropriate QC TSV.
3. Confirm engine/reference/input checksums exist in the manifest.
4. Read `convergence.json`; do not equate engine completion with stability.
5. Open `report.html` for the narrative summary.
6. Use the raw outputs for downstream quantitative work.
7. Preserve the entire output directory, not just the HTML report.

### 14.3 Diagnosing failures

| Symptom | First place to look | Typical cause |
|---|---|---|
| Wrapper exits before engine | terminal message, `run_manifest.json` if present | Invalid flags, unsupported mode/evidence, missing file/reference, failed QC |
| Engine exits nonzero | `pigean_run.log`, `priors_command.txt` | Engine numerical/path/data issue |
| Empty/missing raw output | `pigean_run.log`, manifest execution status | Engine failed before writing the file |
| Many unresolved genes | QC TSV and normalized input | Wrong identifier type, spelling, or gene-map mismatch |
| GWAS missing references | resolved build and `data/` | hg38 requested or reference pack incomplete |
| Docker file not found | mount declarations and container path | Host path not shared/mounted or wrong in-container path |
| Analysis rejected although CLI lists it | compatibility matrix | Adapter exists but end-to-end engine route is blocked/not verified |

### 14.4 Reproducing a run

Use these artifacts together:

- preserved `input/` files;
- `normalized/` files;
- `resolved_config.json`;
- `priors_command.txt`;
- `run_manifest.json` checksums and versions;
- container tag or local Python/NumPy/SciPy versions; and
- genome/gene-set reference checksums.

Because sampling is stochastic, reproduction means the same configuration and
scientifically consistent results within expected variability, not identical
floating-point outputs.

## 15. How to Extend the Pipeline Safely

### 15.1 Adding a new evidence type

1. Trace the exact engine input contract and transformations.
2. Define units, required fields, missing-value behavior, and duplicates.
3. Implement a `BaseEvidenceAdapter` subclass.
4. Add it to the adapter registry and CLI choices.
5. Add configuration/WDL mappings.
6. Add QC, normalization provenance, and fixtures.
7. Declare compatibility in every exposed mode; default to `NOT_VERIFIED`.
8. Run representative end-to-end engine tests.
9. Trace all new output columns before exposing them in reports.
10. Run unit and golden regression suites.

### 15.2 Adding a new model mode

1. Treat the mode as a new estimator, not just a subcommand string.
2. Document required inputs, all defaults, and safe parameter ranges.
3. Trace the engine call chain and output schemas.
4. Implement an `AdvancedMode` subclass.
5. Declare evidence compatibility explicitly.
6. Declare stability applicability and interpretation compatibility.
7. Update parsers if schemas differ.
8. Add mode-specific report language and anti-interpretation.
9. Add failure-manifest coverage for unsupported combinations.
10. Perform end-to-end comparisons without claiming equivalence unless that is
    the actual scientific result.

### 15.3 Adding a genome build

1. Obtain authoritative gene, TSS, and exon location files for the build.
2. Verify their field formats and chromosome conventions against the engine.
3. Add paths under `REFERENCE_PATHS`.
4. Add files to the Docker image and publishing process.
5. Add alias/reference-resolution tests.
6. Run coordinate-based GWAS validation with known loci.
7. Record checksums and provenance.
8. Do not infer or auto-liftover researcher data silently.

### 15.4 Adding or changing interpreted metrics

1. Trace the value from engine calculation through output writing.
2. Record formula, units/scale, direction, stochasticity, and caveats.
3. Classify it in `docs/output_contract.md`.
4. Add it to `INTERPRETATION_CONTRACTS` only if verified.
5. Add parser/report tests and check both text and HTML.
6. Avoid causal or externally calibrated language without evidence.

### 15.5 Release checklist

1. Update `pigean.__version__`.
2. Update `CHANGELOG.md`.
3. Synchronize Docker tags, WDL image tags, and documentation examples.
4. Run 346+ unit tests and the 56-check golden suite.
5. Build and run the container.
6. Verify `engine/priors.py` checksum unless an engine release is intentional.
7. Inspect the full staged diff and exclude generated results/metadata.
8. Commit with the established repository identity.
9. Push code and separately confirm the required container image is published.

## 16. Known Limitations and Technical Debt

### Scientific/runtime limitations

1. **hg19 only in practice.** hg38 aliases exist, but reference files do not.
2. **Z-score mode is blocked.** The frozen engine can hit a SciPy scale-domain
   error for relevant inputs.
3. **Percentile mode is blocked.** The frozen engine references an incorrect
   option attribute.
4. **Naive-priors supports only positive controls.** Other combinations have
   not passed an interpretation-aware end-to-end validation.
5. **No formal MCMC convergence.** Required chain samples are not retained.
6. **Unseeded stochastic output.** Exact repeatability is not available.
7. **Factor/PheWAS/anchors/internal hyperparameters are not exposed.** Their
   contracts are specialized, incomplete, or broken.
8. **No automatic liftover.** Coordinate conversion is the researcher's
   responsibility and must occur before the run.

### Interface/documentation debt

1. The CLI help still shows Z-score as an example and lists both blocked
   adapters as choices. The mode layer rejects them safely, but help text could
   be clearer.
2. `docs/evidence_inputs.md` is the historical Phase 3 document. It predates
   GWAS and the Phase 6 compatibility layer, so its support table is not current.
3. Some early genome-build failures occur before a failure manifest is written.
4. Background use is not explicitly rejected for non-positive-control evidence
   at argument-validation time; the documented and validated use is positive
   controls only.
5. The report caveat may mention increasing `--max-num-iter`, but that setting
   is not exposed as a normal wrapper CLI option in version 6.0.0.
6. WDL exposes parameters corresponding to blocked adapters because it mirrors
   the wider CLI surface; compatibility is still enforced at runtime.

These are maintenance candidates, not permission to weaken the current safety
gates.

## 17. Documentation Authority and Historical Records

Use documents in this order when they disagree:

1. Current executable code and tests.
2. This handover and [How to Run PIGEAN](how_to_run.md).
3. [Advanced Modes](advanced_modes.md), [Output Contract](output_contract.md),
   [Convergence](convergence.md), and [GWAS Input](gwas_input.md).
4. [Architecture](architecture.md) and [Repository Structure](repository_structure.md).
5. Historical records under `docs/phases/`.

Historical documents are valuable for understanding why a decision was made,
but plans may describe targets that were later blocked or revised. In
particular:

- `docs/phases/plan.README.md` is the original roadmap plus implementation
  journal, not the current support contract;
- `docs/phases/execution.md` describes early implementation design;
- `docs/phases/results.README.md` is an early raw-results guide; and
- `docs/evidence_inputs.md` captures Phase 3 before GWAS and Phase 6 mode
  validation.

When adding a feature, update the current contract documents and changelog; do
not rewrite historical records to make history look cleaner.

## 18. Glossary

| Term | Meaning in this project |
|---|---|
| Adapter | Evidence-specific validation, normalization, and engine-argument layer |
| Annotation/gene set | Named collection of genes used as a model feature |
| Background | Eligible positive-control population passed as `--positive-controls-all-in` |
| Candidate gene | Non-input gene prioritized by the model; not a causal claim |
| Combined | Gene prior plus observed evidence, on a log-odds scale |
| `combined_D` | Standard-mode sigmoid of `combined`; model-derived posterior probability |
| Engine | Frozen `engine/priors.py` scientific implementation |
| Evidence type | Researcher input selected by `--analysis` |
| Golden test | Tolerance-based regression specification for the standard baseline |
| HuGE score | Engine-derived gene-level evidence from GWAS/exome association input |
| Manifest | JSON audit record of settings, checksums, validation, execution, and stability |
| Mode | Verified statistical engine pathway selected by `--mode` |
| Normalized input | Preserved engine-ready copy created by the wrapper |
| PIGEAN stability criterion | Engine-specific max-fractional-SEM threshold, not formal MCMC convergence |
| Positive control | Researcher-supplied input/training gene |
| Prior | Gene-set-derived gene log-odds contribution; estimator depends on mode |
| Reference pack | Coordinate and mapping files tied to a genome build |
| Wrapper | CLI and `pigean/` modules surrounding the frozen engine |

---

The most important handover principle is simple: **a feature is not supported
because the engine has a flag or because an adapter can parse its input. It is
supported only when input, engine behavior, output schema, stability
applicability, interpretation, provenance, and regression evidence form one
validated end-to-end contract.**
