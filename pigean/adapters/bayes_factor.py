"""Bayes factor evidence adapter.

Validates, normalizes, and translates a gene-level Bayes factor file
into the --gene-bfs-in argument for priors.py.

Engine contract (traced from priors.py _read_gene_bfs, line 11606):
  - Header required
  - Default gene column: "Gene"
  - Default score column: "log_bf" (natural-log Bayes factor)
  - Whitespace-delimited
  - "NA" values silently skipped by engine
  - Non-numeric values warned and skipped by engine
  - Duplicate genes: last value wins (dict overwrite)
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


class BayesFactorAdapter(BaseEvidenceAdapter):
    """Adapter for gene-level log Bayes factor inputs."""

    ANALYSIS_TYPE = "gene-bayes-factor"
    ENGINE_FLAG = "--gene-bfs-in"
    EVIDENCE_DESCRIPTION = (
        "Gene-level natural-log Bayes factors (header with Gene and log_bf columns)"
    )

    def validate_file(self, input_path, params=None):
        """Validate the BF evidence file.

        Checks file structure, then verifies required columns exist.
        The gene column defaults to "Gene" and the score column to "log_bf"
        unless overridden via params.
        """
        params = params or {}
        gene_col = params.get("gene_column", "Gene")
        score_col = params.get("score_column", "log_bf")

        result = validate_tabular_file(
            input_path,
            required_columns=[gene_col, score_col],
            min_data_rows=1,
        )
        if result.status == ValidationStatus.FAIL:
            return result

        issues = []
        warnings = []
        with open(input_path, "r") as f:
            header = f.readline().strip().split()
            gene_idx = header.index(gene_col)
            score_idx = header.index(score_col)

            row_num = 0
            valid_rows = 0
            for line in f:
                row_num += 1
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
                issues=["No valid data rows with finite numeric log_bf values"],
            )

        if issues:
            return FileValidationResult(
                status=ValidationStatus.PASS_WITH_WARNINGS,
                issues=issues,
            )

        return FileValidationResult(status=ValidationStatus.PASS)

    def normalize(self, input_path, output_dir, params=None):
        """Normalize the BF evidence file.

        1. Preserves original in output_dir/input/
        2. Writes engine-ready file to output_dir/normalized/
           - Keeps only rows with valid gene + finite numeric score
           - Preserves the Gene and log_bf columns the engine expects
           - Reports duplicates (engine uses last-wins, so we preserve all
             and let engine handle it, but warn)
        """
        params = params or {}
        gene_col = params.get("gene_column", "Gene")
        score_col = params.get("score_column", "log_bf")

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, "bayes_factor_original.tsv")
        shutil.copy2(input_path, original_dest)

        warnings = []
        transformations = []

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
                warnings.append(f"Removed row with too few columns")
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

        need_rename = (gene_col != "Gene" or score_col != "log_bf")

        normalized_path = os.path.join(norm_dir, "bayes_factor.tsv")
        with open(normalized_path, "w") as f:
            if need_rename:
                out_header = list(header)
                out_header[gene_idx] = "Gene"
                out_header[score_idx] = "log_bf"
                f.write("\t".join(out_header) + "\n")
                transformations.append(
                    f"Renamed columns: {gene_col}->Gene, {score_col}->log_bf"
                )
            else:
                f.write("\t".join(header) + "\n")

            for line in valid_lines:
                cols = line.strip().split()
                f.write("\t".join(cols) + "\n")

        transformations.append("No scale transformation (values are natural-log BF)")

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
        """Return priors.py CLI arguments for BF input.

        Uses --gene-bfs-in with the engine's default column expectations
        (Gene, log_bf). The normalized file already has these column names.
        """
        return ["--gene-bfs-in", normalized_path]
