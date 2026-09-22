# Repository Inventory — Pre-Reorganization Snapshot

Captured: 2026-09-21, before Phase 5.5 reorganization.

## Root-Level Files

| File | Size | Description |
|------|------|-------------|
| `priors.py` | 952 KB | Frozen PIGEAN engine (16,567 lines) |
| `run_pigean.py` | 22 KB | CLI entry point |
| `plan.README.md` | 131 KB | Implementation plan (Phases 0–7 design) |
| `execution.md` | 62 KB | Implementation journal |
| `results.README.md` | 14 KB | Results interpretation guide |
| `README.md` | 3.5 KB | Project readme |
| `LICENSE` | 1.5 KB | BSD 3-Clause (Broad Institute) |
| `Dockerfile` | 570 B | Docker build file |
| `gclouddock.sh` | 453 B | GCR push script (stale tag 1.0.0) |
| `.gitignore` | 6 B | Only `.idea/` |

## Root-Level Stray Output Files

These are from ad-hoc engine runs during development. Not referenced by any code.

| File | Size | Date | SHA-256 |
|------|------|------|---------|
| `gs.out` | 2.9 MB | Sep 11 | `843ee227…` |
| `gss.out` | 6.8 MB | Sep 11 | `572248…` |
| `ggss.out` | 19.9 MB | Sep 11 | `26ee40…` |
| `p.out` | 936 B | Sep 11 | `051048…` |
| `gs.full.out` | 2.9 MB | Sep 9 | `a100b0…` |
| `gss.full.out` | 6.8 MB | Sep 9 | `acfba6…` |
| `ggss.full.out` | 20.1 MB | Sep 9 | `c5e5e1…` |
| `p.full.out` | 936 B | Sep 9 | `e4e92e…` |

## Directories

| Directory | Contents |
|-----------|----------|
| `pigean/` | Wrapper Python package (10 modules + 8 adapters) |
| `data/` | 6 hg19 reference/annotation files (25 MB) |
| `tests/` | 9 test files, golden test, 3 fixture directories |
| `wdl/` | 2 WDL workflow files |
| `ex/` | Example gene list (1 file) |
| `raw/` | Original Broad cluster command (1 file) |
| `outputs/` | Phase 0–5 development artifacts (85 MB) |
| `runs/` | Reference wrapper run (8.7 MB) |
| `docs/` | 4 documentation files |
| `__pycache__/` | Python bytecode cache |
| `.pytest_cache/` | pytest cache |

## Frozen Engine Checksum

```
da66b0c75f4554d7256daa5a2078e05113af7d746e5102d3b232ccbece2fb165  priors.py
```

## Reference Data Checksums

```
bf5635ecd4f485b03d6f266cdc61d37a19e22c3a1afe272323e1889f14c54597  data/gene_set_list_mouse_2024.txt
42ed39dc4ff8ad51fab3c031a3fdfa88709abdc0fa7784dd58410011169d0f4d  data/gene_set_list_msigdb_nohp.txt
abc7bbce31c4b2fe05ff2d410d9ca1081538d6852d65c11376449b5a9f6aaeed  data/NCBI37.3.plink.gene.exons.loc
0f37c84a9ba460ac88c416863146c8e05b309d35cadf7ccff596be610365242b  data/NCBI37.3.plink.gene.loc
65a9fd75392b519363491e4a26c99a2fbcffc9501132057ea5481f2ebcd49132  data/portal_gencode.gene.map
1e6aa427d730f0cb75bc3b4df5a1e44ed1eb8b02a8b6f21f387b4ff8be3a0127  data/refGene_hg19_TSS.subset.loc
```

## Test Baseline

- Unit tests: 323 PASS, 0 FAIL
- Golden regression: 56 PASS, 0 FAIL
