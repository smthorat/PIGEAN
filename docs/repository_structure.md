# Repository Structure

## Directory Layout

```
rock-pigean/
├── run_pigean.py              # CLI entry point (wrapper)
├── pigean/                    # Wrapper Python package
│   ├── __init__.py            #   Package version (5.0.0)
│   ├── config.py              #   Configuration resolution
│   ├── engine.py              #   Engine command building + execution
│   ├── validation.py          #   Input file + gene QC
│   ├── references.py          #   Path resolution (Docker / local)
│   ├── parsers.py             #   Output file parsing (gs/gss/ggss/p.out)
│   ├── convergence.py         #   PIGEAN stability assessment
│   ├── interpretation.py      #   Results interpretation
│   ├── report.py              #   Text + HTML report generation
│   ├── manifest.py            #   Run provenance manifest
│   └── adapters/              #   Evidence type adapters
│       ├── __init__.py        #     Registry + base types
│       ├── base.py            #     BaseEvidenceAdapter
│       ├── positive_controls.py
│       ├── bayes_factor.py
│       ├── zscore.py
│       ├── percentile.py
│       ├── exome.py
│       └── gwas.py
├── engine/
│   └── priors.py              # Frozen PIGEAN engine (DO NOT MODIFY)
├── data/                      # Scientific reference data (hg19)
│   ├── gene_set_list_mouse_2024.txt    # Mouse phenotype gene sets
│   ├── gene_set_list_msigdb_nohp.txt   # MSigDB gene sets
│   ├── portal_gencode.gene.map         # Gene name mapping
│   ├── NCBI37.3.plink.gene.loc         # Gene locations
│   ├── NCBI37.3.plink.gene.exons.loc   # Exon-level locations
│   └── refGene_hg19_TSS.subset.loc     # TSS locations
├── tests/
│   ├── test_config.py
│   ├── test_validation.py
│   ├── test_engine_command.py
│   ├── test_gene_sets.py
│   ├── test_evidence_adapters.py
│   ├── test_gwas_adapter.py
│   ├── test_parsers.py
│   ├── test_convergence.py
│   ├── test_interpretation.py
│   ├── golden/                # Golden regression test
│   │   ├── run_golden_test.sh
│   │   ├── checkpoints.json
│   │   ├── input/
│   │   └── expected/
│   └── fixtures/              # Test fixtures
│       ├── evidence/
│       ├── gwas/
│       └── outputs/
├── workflows/                 # WDL for Terra/Cromwell
│   ├── rock_pigean.wdl
│   └── grep_genes.wdl
├── examples/
│   └── gene_list              # Example positive-controls input
├── scripts/
│   ├── gclouddock.sh          # Docker build + GCR push
│   └── legacy/
│       └── pigean.sh          # Original Broad cluster command
├── results/
│   ├── development/           # Phase validation artifacts
│   │   ├── phase0_golden/     # Stochastic tolerance estimation
│   │   ├── phase1_wrapper/    # Wrapper equivalence tests
│   │   ├── phase2_custom_gene_sets/
│   │   ├── phase3_evidence/
│   │   ├── phase4_gwas/
│   │   ├── phase5_interpretation/
│   │   ├── phase5_report/
│   │   └── scratch_outputs/   # Legacy ad-hoc engine output files
│   └── runs/
│       └── phase1_run/        # Reference complete wrapper run
├── docs/
│   ├── architecture.md        # System architecture
│   ├── repository_structure.md  # This file
│   ├── convergence.md         # Stability assessment docs
│   ├── evidence_inputs.md     # Evidence input formats
│   ├── gwas_input.md          # GWAS input format
│   ├── output_contract.md     # Output file column contracts
│   ├── path_dependency_audit.md  # Path dependency analysis
│   ├── repository_inventory_before.md  # Pre-reorganization snapshot
│   └── phases/                # Phase history documents
│       ├── plan.README.md     # Original implementation plan
│       ├── execution.md       # Implementation journal
│       └── results.README.md  # Results interpretation guide
├── Dockerfile
├── README.md
├── LICENSE                    # BSD 3-Clause (Broad Institute)
├── CHANGELOG.md
└── .gitignore
```

## Design Decisions

### `data/` kept at root (not renamed to `resources/`)

The name `data/` is standard and clear. Renaming would require updating ~28
path references across `pigean/references.py`, `Dockerfile`,
`tests/golden/run_golden_test.sh`, and `tests/test_engine_command.py` with
minimal clarity gain. Deferred.

### `run_pigean.py` and `pigean/` kept at root (not moved to `code/`)

Moving the Python package would break every import statement. Moving the entry
point would require Docker ENTRYPOINT, WDL command, and `sys.path` changes.
Both are already in standard Python project positions.

### `tests/` kept at root

Standard Python project convention. All test files use
`sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))` to import
from `pigean/`, which works because `tests/` is one level below repo root.
