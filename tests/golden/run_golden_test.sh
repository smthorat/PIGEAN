#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_golden_test.sh — Phase 0 Golden Regression Test
#
# Runs PIGEAN with the frozen test input and validates the output against
# the golden expected values using empirically derived tolerances.
#
# Usage:
#   bash tests/golden/run_golden_test.sh [--keep]              # run engine fresh
#   bash tests/golden/run_golden_test.sh <output_dir>          # validate existing outputs
#   bash tests/golden/run_golden_test.sh <output_dir> --keep   # (--keep is no-op here)
#
# Outputs:
#   - Prints PASS/FAIL for each check
#   - Exit 0 if all checks pass, 1 if any fail
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GOLDEN_DIR="$SCRIPT_DIR"
EXPECTED_DIR="$GOLDEN_DIR/expected"
CHECKPOINTS="$GOLDEN_DIR/checkpoints.json"

# Parse args
KEEP=false
EXISTING_DIR=""
for arg in "$@"; do
    case $arg in
        --keep) KEEP=true ;;
        *)
            if [ -d "$arg" ]; then
                EXISTING_DIR="$arg"
            fi
            ;;
    esac
done

if [ -n "$EXISTING_DIR" ]; then
    TMPDIR="$EXISTING_DIR"
    KEEP=true  # never delete a user-supplied directory
    SKIP_ENGINE=true
else
    TMPDIR=$(mktemp -d /tmp/pigean_golden_test.XXXXXX)
    trap '[ "$KEEP" = false ] && rm -rf "$TMPDIR"' EXIT
    SKIP_ENGINE=false
fi

RUN_LABEL="${EXISTING_DIR:-fresh run}"
echo "════════════════════════════════════════════════════════════════"
echo "  PIGEAN Golden Regression Test"
echo "  $(date)"
echo "  Source: $RUN_LABEL"
echo "  Output: $TMPDIR"
echo "════════════════════════════════════════════════════════════════"

PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

pass() { echo "  ✓ PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  ✗ FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }
warn() { echo "  ⚠ WARN: $1"; WARN_COUNT=$((WARN_COUNT + 1)); }

# ── Verify prerequisites ───────────────────────────────────────────────────
echo ""
echo "Checking prerequisites..."

if [ ! -f "$CHECKPOINTS" ]; then
    echo "ERROR: checkpoints.json not found at $CHECKPOINTS"
    exit 1
fi

if [ ! -f "$GOLDEN_DIR/input/positive_controls.txt" ]; then
    echo "ERROR: Golden input not found"
    exit 1
fi

for gfile in gs.golden.out gss.golden.out ggss.golden.out p.golden.out; do
    if [ ! -f "$EXPECTED_DIR/$gfile" ]; then
        echo "ERROR: Expected output $gfile not found"
        exit 1
    fi
done

# ── Run PIGEAN engine (or skip if validating existing outputs) ──────────────
echo ""

if [ "$SKIP_ENGINE" = true ]; then
    echo "Validating existing outputs in: $TMPDIR"
    # Verify required files exist
    for f in gs.out gss.out p.out; do
        if [ ! -f "$TMPDIR/$f" ]; then
            echo "ERROR: $TMPDIR/$f not found"
            exit 1
        fi
    done
    pass "Using pre-computed outputs"
else
    echo "Running PIGEAN engine..."

    BASEDIR="$REPO_ROOT"
    python3 -u "$BASEDIR/engine/priors.py" gibbs \
        --X-in "$BASEDIR/data/gene_set_list_mouse_2024.txt" \
        --X-in "$BASEDIR/data/gene_set_list_msigdb_nohp.txt" \
        --gene-map-in "$BASEDIR/data/portal_gencode.gene.map" \
        --max-num-gene-sets 5000 \
        --gene-stats-out "$TMPDIR/gs.out" \
        --gene-set-stats-out "$TMPDIR/gss.out" \
        --gene-gene-set-stats "$TMPDIR/ggss.out" \
        --params-out "$TMPDIR/p.out" \
        --debug-level 3 \
        --positive-controls-in "$GOLDEN_DIR/input/positive_controls.txt" \
        --gene-loc-file "$BASEDIR/data/NCBI37.3.plink.gene.loc" \
        --gene-loc-file-huge "$BASEDIR/data/refGene_hg19_TSS.subset.loc" \
        --exons-loc-file-huge "$BASEDIR/data/NCBI37.3.plink.gene.exons.loc" \
        --gene-filter-value 1 --gene-set-filter-value 0.01 \
        > "$TMPDIR/pigean_run.log" 2>&1

    ENGINE_EXIT=$?

    if [ $ENGINE_EXIT -ne 0 ]; then
        fail "Engine exited with code $ENGINE_EXIT"
        echo "  See log: $TMPDIR/pigean_run.log"
        echo ""
        echo "RESULT: FAIL ($FAIL_COUNT failures)"
        exit 1
    fi
    pass "Engine completed successfully (exit code 0)"
fi

# ── Run validation via Python ──────────────────────────────────────────────
echo ""
echo "Validating outputs..."

python3 - "$TMPDIR" "$CHECKPOINTS" "$EXPECTED_DIR" << 'PYEOF'
import sys
import os
import json
import csv

tmpdir = sys.argv[1]
checkpoints_path = sys.argv[2]
expected_dir = sys.argv[3]

with open(checkpoints_path) as f:
    ckpt = json.load(f)

passes = []
fails = []
warns = []

def check_pass(msg):
    passes.append(msg)
    print(f"  ✓ PASS: {msg}")

def check_fail(msg):
    fails.append(msg)
    print(f"  ✗ FAIL: {msg}")

def check_warn(msg):
    warns.append(msg)
    print(f"  ⚠ WARN: {msg}")

# ── Structural checks ──────────────────────────────────────────────────

# Row counts
for fname, key in [("gs.out", "gs_data_rows"), ("gss.out", "gss_data_rows"), ("p.out", "p_data_rows")]:
    fpath = os.path.join(tmpdir, fname)
    if not os.path.exists(fpath):
        check_fail(f"{fname} not produced")
        continue
    with open(fpath) as f:
        data_rows = sum(1 for _ in f) - 1  # exclude header
    expected = ckpt["structural"][key]
    max_key = key + "_max"
    expected_max = ckpt["structural"].get(max_key, expected)

    if expected <= data_rows <= expected_max:
        check_pass(f"{fname} row count: {data_rows} (expected {expected}–{expected_max})")
    elif data_rows == expected:
        check_pass(f"{fname} row count: {data_rows} (expected {expected})")
    else:
        pct_diff = abs(data_rows - expected) / expected * 100
        if pct_diff < 1:
            check_warn(f"{fname} row count: {data_rows} (expected {expected}, {pct_diff:.2f}% off)")
        else:
            check_fail(f"{fname} row count: {data_rows} (expected {expected}, {pct_diff:.2f}% off)")

# Header check
for fname, key in [("gs.out", "gs_header"), ("gss.out", "gss_header")]:
    fpath = os.path.join(tmpdir, fname)
    if not os.path.exists(fpath):
        continue
    with open(fpath) as f:
        header = f.readline().strip().split("\t")
    expected_header = ckpt["structural"][key]
    if header == expected_header:
        check_pass(f"{fname} header matches expected")
    else:
        check_fail(f"{fname} header mismatch: got {header[:5]}..., expected {expected_header[:5]}...")

# ── Gene priors ────────────────────────────────────────────────────────

gs_path = os.path.join(tmpdir, "gs.out")
if os.path.exists(gs_path):
    with open(gs_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        gs_data = {r["Gene"]: r for r in reader}

    for gene, spec in ckpt.get("gene_priors", {}).items():
        if gene not in gs_data:
            check_fail(f"Gene {gene} not found in gs.out")
            continue

        actual = float(gs_data[gene]["prior"])
        expected = spec["expected"]
        tolerance = spec["tolerance"]

        if abs(actual - expected) <= tolerance:
            check_pass(f"{gene} prior={actual:.6f} (expected {expected:.6f} ± {tolerance})")
        else:
            check_fail(f"{gene} prior={actual:.6f} (expected {expected:.6f} ± {tolerance}, off by {abs(actual-expected):.6f})")

    # ── Rank correlation with golden ───────────────────────────────────
    # Load golden gs.out
    golden_gs_path = os.path.join(expected_dir, "gs.golden.out")
    if os.path.exists(golden_gs_path):
        with open(golden_gs_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            golden_gs = {r["Gene"]: r for r in reader}

        # Compute Spearman for common genes
        common = sorted(set(gs_data.keys()) & set(golden_gs.keys()))
        if len(common) > 100:
            from scipy.stats import spearmanr
            for col in ["prior", "combined", "log_bf"]:
                try:
                    actual_vals = [float(gs_data[g][col]) for g in common]
                    golden_vals = [float(golden_gs[g][col]) for g in common]
                    rho, _ = spearmanr(actual_vals, golden_vals)

                    min_rho = ckpt.get("rank_stability", {}).get(col, {}).get("min_spearman_rho", 0.95)
                    if rho >= min_rho:
                        check_pass(f"Spearman rho ({col}): {rho:.6f} ≥ {min_rho}")
                    else:
                        check_fail(f"Spearman rho ({col}): {rho:.6f} < {min_rho}")
                except Exception as e:
                    check_warn(f"Could not compute Spearman for {col}: {e}")

    # ── Top gene overlap ──────────────────────────────────────────────
    golden_top = ckpt.get("golden_top_non_input_genes", [])
    if golden_top:
        input_genes = {"RMI2", "LITAF", "SLCO1B1", "BRD2", "SOX2", "LEP", "PLINK1", "ADIPOQ", "LECT2", "CD300LG"}
        non_input = [(g, float(d["prior"])) for g, d in gs_data.items() if g not in input_genes]
        non_input.sort(key=lambda x: x[1], reverse=True)
        actual_top20 = set(g for g, _ in non_input[:20])
        golden_top20 = set(golden_top[:20])
        overlap = actual_top20 & golden_top20

        min_overlap = ckpt.get("top_gene_overlap", {}).get("top_20", {}).get("min_overlap", 10)
        if len(overlap) >= min_overlap:
            check_pass(f"Top-20 non-input gene overlap: {len(overlap)}/20 ≥ {min_overlap}")
        else:
            check_fail(f"Top-20 non-input gene overlap: {len(overlap)}/20 < {min_overlap}")

# ── Gene set / pathway checks ─────────────────────────────────────────

gss_path = os.path.join(tmpdir, "gss.out")
golden_gss_path = os.path.join(expected_dir, "gss.golden.out")
if os.path.exists(gss_path) and os.path.exists(golden_gss_path):
    with open(gss_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        gss_data = {r["Gene_Set"]: r for r in reader}

    with open(golden_gss_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        golden_gss = {r["Gene_Set"]: r for r in reader}

    # Spearman for gene set betas
    common_gs = sorted(set(gss_data.keys()) & set(golden_gss.keys()))
    if len(common_gs) > 10:
        from scipy.stats import spearmanr
        try:
            actual_betas = [float(gss_data[gs]["beta"]) for gs in common_gs]
            golden_betas = [float(golden_gss[gs]["beta"]) for gs in common_gs]
            rho, _ = spearmanr(actual_betas, golden_betas)

            min_rho = ckpt.get("rank_stability", {}).get("gene_set_beta", {}).get("min_spearman_rho", 0.95)
            if rho >= min_rho:
                check_pass(f"Gene set beta Spearman rho: {rho:.6f} ≥ {min_rho}")
            else:
                check_fail(f"Gene set beta Spearman rho: {rho:.6f} < {min_rho}")
        except Exception as e:
            check_warn(f"Could not compute gene set Spearman: {e}")

    # Top pathway overlap
    golden_top_gs = ckpt.get("golden_top_gene_sets", [])
    if golden_top_gs:
        actual_sorted = sorted(gss_data.items(), key=lambda x: float(x[1]["beta"]), reverse=True)
        actual_top20 = set(gs for gs, _ in actual_sorted[:20])
        golden_top20 = set(golden_top_gs[:20])
        overlap = actual_top20 & golden_top20

        min_overlap = ckpt.get("top_pathway_overlap", {}).get("top_20", {}).get("min_overlap", 10)
        if len(overlap) >= min_overlap:
            check_pass(f"Top-20 pathway overlap: {len(overlap)}/20 ≥ {min_overlap}")
        else:
            check_fail(f"Top-20 pathway overlap: {len(overlap)}/20 < {min_overlap}")

# ── Model parameters ──────────────────────────────────────────────────

p_path = os.path.join(tmpdir, "p.out")
if os.path.exists(p_path):
    with open(p_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        p_data = {}
        for r in reader:
            p_data[(r["Parameter"], r["Version"])] = r["Value"]

    for key, spec in ckpt.get("model_parameters", {}).items():
        # key format: "param_name_vN"
        parts = key.rsplit("_v", 1)
        if len(parts) != 2:
            continue
        param, version = parts[0], parts[1]
        actual_val = p_data.get((param, version))

        if actual_val is None:
            check_warn(f"Parameter {param} v{version} not found in p.out")
            continue

        if spec.get("exact", False):
            if actual_val == spec["expected"]:
                check_pass(f"Parameter {param} v{version} = {actual_val} (exact match)")
            else:
                check_fail(f"Parameter {param} v{version} = {actual_val} (expected {spec['expected']})")
        else:
            try:
                actual_num = float(actual_val)
                expected_num = float(spec["expected"])
                tolerance = float(spec["tolerance"])
                if abs(actual_num - expected_num) <= tolerance:
                    check_pass(f"Parameter {param} v{version} = {actual_num:.6g} (within ±{tolerance})")
                else:
                    check_fail(f"Parameter {param} v{version} = {actual_num:.6g} (expected {expected_num:.6g} ±{tolerance})")
            except ValueError:
                if actual_val == spec.get("expected"):
                    check_pass(f"Parameter {param} v{version} = {actual_val}")
                else:
                    check_fail(f"Parameter {param} v{version} = {actual_val} (expected {spec.get('expected')})")

# ── Summary ───────────────────────────────────────────────────────────
print()
print("═" * 60)
total = len(passes) + len(fails) + len(warns)
print(f"  Results: {len(passes)} PASS, {len(fails)} FAIL, {len(warns)} WARN (total {total})")
if fails:
    print(f"  OVERALL: FAIL")
    sys.exit(1)
else:
    print(f"  OVERALL: PASS")
    sys.exit(0)
PYEOF

VALIDATION_EXIT=$?

echo ""
echo "════════════════════════════════════════════════════════════════"
if [ $VALIDATION_EXIT -eq 0 ]; then
    echo "  GOLDEN TEST: PASS"
    echo "════════════════════════════════════════════════════════════════"
    if [ "$KEEP" = true ]; then
        echo "  Output kept at: $TMPDIR"
    fi
    exit 0
else
    echo "  GOLDEN TEST: FAIL"
    echo "════════════════════════════════════════════════════════════════"
    echo "  Output at: $TMPDIR"
    KEEP=true  # keep output on failure for debugging
    exit 1
fi
