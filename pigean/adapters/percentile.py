"""Gene percentile evidence adapter.

Validates, normalizes, and translates a gene-level percentile/ranking file
into the --gene-percentiles-in argument for priors.py.

Engine contract (traced from priors.py _read_gene_percentiles, line 11843):
  - Header required
  - Gene ID column MUST be specified (--gene-percentiles-id-col)
  - Value column MUST be specified (--gene-percentiles-value-col)
  - Default: lower value = better rank (gene_percentiles_higher_is_better=False)
  - Conversion: rank → quantile → inverse normal CDF → log-BF
  - Whitespace-delimited
  - "NA" values silently skipped
  - Non-numeric values warned and skipped
  - Duplicate genes: last value wins
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


class PercentileAdapter(BaseEvidenceAdapter):
    """Adapter for gene-level percentile/ranking inputs."""

    ANALYSIS_TYPE = "gene-percentile"
    ENGINE_FLAG = "--gene-percentiles-in"
    EVIDENCE_DESCRIPTION = (
        "Gene-level ranking scores. Engine converts ranks to log-BFs "
        "via inverse normal transformation. Requires --gene-column and "
        "--score-column."
    )

    def validate_file(self, input_path, params=None):
        """Validate the percentile evidence file.

        Engine requires both gene ID and value column names explicitly.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        if not gene_col:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["--gene-column is required for gene-percentile analysis"],
            )
        if not score_col:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["--score-column is required for gene-percentile analysis"],
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
        """Normalize the percentile evidence file.

        Preserves original, writes engine-ready file with valid rows only.
        No value transformation — the engine handles the full
        rank → quantile → inverse normal → log-BF conversion internally.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, "percentile_original.tsv")
        shutil.copy2(input_path, original_dest)

        warnings = []
        transformations = [
            "No wrapper-level transformation. Engine converts: "
            "rank → quantile → inverse_normal(quantile) → log-BF"
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

        normalized_path = os.path.join(norm_dir, "percentile.tsv")
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
        """Return priors.py CLI arguments for percentile input.

        Engine REQUIRES --gene-percentiles-id-col and --gene-percentiles-value-col.
        Optionally supports --gene-percentiles-higher-is-better flag.
        """
        params = params or {}
        gene_col = params.get("gene_column")
        score_col = params.get("score_column")

        args = [
            "--gene-percentiles-in", normalized_path,
            "--gene-percentiles-id-col", gene_col,
            "--gene-percentiles-value-col", score_col,
        ]

        if params.get("higher_is_better", False):
            args.append("--gene-percentiles-higher-is-better")

        return args
