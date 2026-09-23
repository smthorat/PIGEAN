"""
pigean/parsers.py — Parse PIGEAN engine output files.

Phase 5: Output parsers for gs.out, gss.out, ggss.out, p.out.
stdlib only — no pandas or external dependencies.
"""

import os
from collections import defaultdict


# Columns that should stay as strings (never converted to numeric).
_STRING_COLUMNS = {"Gene", "Gene_Set", "label", "Chrom", "gene_set", "batch"}

# Columns that should be converted to int (not float).
_INT_COLUMNS = {"N", "Start", "End", "Version"}


def _safe_convert(value, column_name):
    """Convert a string value to the appropriate Python type.

    Returns str for string columns, int for integer columns,
    float for everything else. Falls back to raw string on failure.
    """
    if column_name in _STRING_COLUMNS:
        return value
    if column_name in _INT_COLUMNS:
        try:
            return int(value)
        except (ValueError, TypeError):
            try:
                return int(float(value))
            except (ValueError, TypeError):
                return value
    # Default: try float
    try:
        return float(value)
    except (ValueError, TypeError):
        return value


def _parse_tsv(path):
    """Read a tab-separated file, returning (columns, rows).

    Returns None if the file doesn't exist or is empty.
    Each row is a list of raw string values (no conversion yet).
    """
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        header_line = f.readline()
        if not header_line.strip():
            return None
        columns = header_line.rstrip("\n").split("\t")
        rows = []
        for line in f:
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            rows.append(stripped.split("\t"))
    return columns, rows


def _rows_to_dicts(columns, rows):
    """Convert raw string rows into list of dicts with type conversion."""
    result = []
    ncols = len(columns)
    for row in rows:
        d = {}
        for i, col in enumerate(columns):
            raw = row[i] if i < len(row) else ""
            d[col] = _safe_convert(raw, col)
        # Store any extra tab-separated fields beyond the header.
        for i in range(ncols, len(row)):
            d[f"_extra_{i}"] = row[i]
        result.append(d)
    return result


# ── gs.out ─────────────────────────────────────────────────────────────

def parse_gene_stats(gs_path):
    """Parse gs.out into a structured dict.

    Returns None if the file doesn't exist or is empty.
    Genes are sorted by combined_D descending (highest first).
    """
    parsed = _parse_tsv(gs_path)
    if parsed is None:
        return None

    columns, rows = parsed
    genes = _rows_to_dicts(columns, rows)

    has_positive_control = "positive_control" in columns
    huge_cols = {"huge_score", "huge_score_gwas", "huge_score_exomes",
                 "huge_score_gwas_uncorrected"}
    has_huge_scores = bool(huge_cols & set(columns))

    # Standard mode sorts by combined_D. Non-outer-Gibbs modes do not emit
    # combined_D, so fall back to the verified combined log-odds and then prior.
    def sort_key(g):
        for column in ("combined_D", "combined", "prior"):
            v = g.get(column)
            if isinstance(v, (int, float)):
                return -v
        return float("inf")

    genes.sort(key=sort_key)

    return {
        "columns": columns,
        "row_count": len(genes),
        "has_positive_control": has_positive_control,
        "has_huge_scores": has_huge_scores,
        "genes": genes,
    }


# ── gss.out ────────────────────────────────────────────────────────────

def parse_gene_set_stats(gss_path):
    """Parse gss.out into a structured dict.

    Returns None if the file doesn't exist or is empty.
    Gene sets are sorted by |beta| descending.
    """
    parsed = _parse_tsv(gss_path)
    if parsed is None:
        return None

    columns, rows = parsed
    gene_sets = _rows_to_dicts(columns, rows)

    def sort_key(gs):
        v = gs.get("beta")
        if isinstance(v, (int, float)):
            return -abs(v)
        return float("inf")

    gene_sets.sort(key=sort_key)

    return {
        "columns": columns,
        "row_count": len(gene_sets),
        "gene_sets": gene_sets,
    }


# ── ggss.out ───────────────────────────────────────────────────────────

def parse_gene_gene_set_stats(ggss_path):
    """Parse ggss.out into a structured dict with bidirectional indexes.

    Returns None if the file doesn't exist or is empty.
    Entries are stored as-is; by_gene and by_gene_set provide fast lookup.
    """
    parsed = _parse_tsv(ggss_path)
    if parsed is None:
        return None

    columns, rows = parsed
    entries = _rows_to_dicts(columns, rows)

    by_gene = defaultdict(list)
    by_gene_set = defaultdict(list)

    for entry in entries:
        gene = entry.get("Gene", "")
        gs = entry.get("gene_set", "")
        if gene:
            by_gene[gene].append(entry)
        if gs:
            by_gene_set[gs].append(entry)

    return {
        "columns": columns,
        "row_count": len(entries),
        "entries": entries,
        "by_gene": dict(by_gene),
        "by_gene_set": dict(by_gene_set),
    }


# ── p.out ──────────────────────────────────────────────────────────────

def parse_params(p_path):
    """Parse p.out (Parameter / Version / Value TSV).

    Multi-version params (same name, different versions) are stored as
    lists. Single-version params are stored as scalars.

    Returns None if the file doesn't exist or is empty.
    """
    parsed = _parse_tsv(p_path)
    if parsed is None:
        return None

    columns, rows = parsed
    raw_rows = []
    # Collect all (param, version, value) triples.
    versions = defaultdict(list)  # param_name -> [(version_int, value)]

    for row in rows:
        param = row[0] if len(row) > 0 else ""
        ver_str = row[1] if len(row) > 1 else "1"
        val_str = row[2] if len(row) > 2 else ""
        raw_rows.append((param, ver_str, val_str))

        # Convert version to int.
        try:
            ver = int(ver_str)
        except (ValueError, TypeError):
            ver = 1

        # Convert value — try int first (for counts), then float, then str.
        value = val_str
        if val_str in ("True", "False"):
            value = val_str == "True"
        elif val_str == "None":
            value = None
        else:
            try:
                # Prefer int for whole numbers.
                fval = float(val_str)
                if fval == int(fval) and "." not in val_str and "e" not in val_str.lower():
                    value = int(fval)
                else:
                    value = fval
            except (ValueError, TypeError):
                value = val_str

        versions[param].append((ver, value))

    # Build params dict: single-version → scalar, multi-version → list.
    params = {}
    for param, ver_vals in versions.items():
        ver_vals.sort(key=lambda x: x[0])
        if len(ver_vals) == 1:
            params[param] = ver_vals[0][1]
        else:
            params[param] = [v for _, v in ver_vals]

    return {
        "params": params,
        "raw_rows": raw_rows,
    }
