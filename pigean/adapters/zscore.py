"""Gene Z-score evidence adapter.

Validates, normalizes, and translates a gene-level Z-score file
into the --gene-zs-in argument for priors.py.

Engine contract (traced from priors.py _read_gene_zs, line 11730):
  - Header required
  - Gene ID column MUST be specified (--gene-zs-id-col)
  - Value column MUST be specified (--gene-zs-value-col)
  - Default mode (use_zs_as_log_odds=True): treats values as raw log-odds,
    NOT statistical Z-scores. Conversion: bf = z - mean(z) + background_log_bf
  - Whitespace-delimited
  - "NA" values silently skipped
  - Non-numeric values warned and skipped
  - Duplicate genes: last value wins

IMPORTANT SEMANTIC NOTE: The engine's "Z-score" mode treats input values
as raw log-odds scores, mean-centers them, and shifts by background_log_bf.
This is NOT the same as statistical Z-scores (beta/SE). The wrapper documents
this clearly and uses the term "gene-scores" to avoid confusion.
"""

import os
import shutil
import math

from pigean.adapters import FileValidationResult, ValidationStatus
from pigean.adapters.base import (
    BaseEvidenceAdapter,
    EvidenceNormalizationResult,
    validate_tabular_file,
)


class ZScoreAdapter(BaseEvidenceAdapter):
    """Adapter for gene-level Z-score/log-odds inputs."""

    ANALYSIS_TYPE = "gene-z-score"
    ENGINE_FLAG = "--gene-zs-in"
    EVIDENCE_DESCRIPTION = (
        "Gene-level scores treated as log-odds by the engine. "
        "Values are mean-centered then shifted by the background log-BF. "
        "Requires --gene-column and --score-column."
    )

    def validate_file(self, input_path, params=None):
        """Validate the Z-score evidence file.

        Engine requires both gene ID and value column names explicitly.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        if not gene_col:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["--gene-column is required for gene-z-score analysis"],
            )
        if not score_col:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["--score-column is required for gene-z-score analysis"],
            )

        result = validate_tabular_file(
            input_path,
            required_columns=[gene_col, score_col],
            min_data_rows=1,
        )
        if result.status == ValidationStatus.FAIL:
            return result

        issues = []
        with open(input_path, "r") as f:
            header = f.readline().strip().split()
            gene_idx = header.index(gene_col)
            score_idx = header.index(score_col)

            valid_rows = 0
            for row_num, line in enumerate(f, 1):
                cols = line.strip().split()
                if not cols:
                    continue
                if score_idx >= len(cols) or gene_idx >= len(cols):
                    issues.append(f"Row {row_num}: too few columns")
                    continue
                val = cols[score_idx]
                if val == "NA":
                    continue
                try:
                    v = float(val)
                    if math.isnan(v) or math.isinf(v):
                        issues.append(f"Row {row_num}: non-finite value {val}")
                        continue
                except ValueError:
                    issues.append(f"Row {row_num}: non-numeric score '{val}'")
                    continue
                valid_rows += 1

        if valid_rows == 0:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["No valid data rows with finite numeric score values"],
            )

        if issues:
            return FileValidationResult(
                status=ValidationStatus.PASS_WITH_WARNINGS,
                issues=issues,
            )

        return FileValidationResult(status=ValidationStatus.PASS)

    def normalize(self, input_path, output_dir, params=None):
        """Normalize the Z-score evidence file.

        Preserves original, writes engine-ready file with valid rows only.
        No value transformation — the engine handles mean-centering internally.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, "zscore_original.tsv")
        shutil.copy2(input_path, original_dest)

        warnings = []
        transformations = [
            "No wrapper-level transformation. Engine treats values as log-odds: "
            "bf = value - mean(values) + background_log_bf"
        ]

        with open(input_path, "r") as f:
            header_line = f.readline()
            header = header_line.strip().split()
            gene_idx = header.index(gene_col)
            score_idx = header.index(score_col)
            data_lines = f.readlines()

        original_rows = len(data_lines)
        removed_rows = 0
        blanks = 0
        duplicates = 0
        seen_genes = set()
        valid_lines = []

        for line in data_lines:
            cols = line.strip().split()
            if not cols:
                blanks += 1
                removed_rows += 1
                continue
            if gene_idx >= len(cols) or score_idx >= len(cols):
                removed_rows += 1
                continue
            gene = cols[gene_idx]
            val = cols[score_idx]
            if val == "NA":
                removed_rows += 1
                continue
            try:
                v = float(val)
                if math.isnan(v) or math.isinf(v):
                    removed_rows += 1
                    warnings.append(f"Removed {gene}: non-finite value {val}")
                    continue
            except ValueError:
                removed_rows += 1
                warnings.append(f"Removed {gene}: non-numeric score '{val}'")
                continue

            if gene in seen_genes:
                duplicates += 1
                warnings.append(f"Duplicate gene: {gene} (engine uses last occurrence)")
            seen_genes.add(gene)
            valid_lines.append(line)

        normalized_path = os.path.join(norm_dir, "zscore.tsv")
        with open(normalized_path, "w") as f:
            f.write("\t".join(header) + "\n")
            for line in valid_lines:
                cols = line.strip().split()
                f.write("\t".join(cols) + "\n")

        return EvidenceNormalizationResult(
            original_path=original_dest,
            normalized_path=normalized_path,
            original_rows=original_rows,
            normalized_rows=len(valid_lines),
            removed_rows=removed_rows,
            duplicates_removed=duplicates,
            blanks_removed=blanks,
            transformations=transformations,
            warnings=warnings,
        )

    def build_engine_args(self, normalized_path, params=None):
        """Return priors.py CLI arguments for Z-score input.

        Engine REQUIRES --gene-zs-id-col and --gene-zs-value-col.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        args = [
            "--gene-zs-in", normalized_path,
            "--gene-zs-id-col", gene_col,
            "--gene-zs-value-col", score_col,
        ]
        return args
