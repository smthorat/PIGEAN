"""Exome/rare-variant association evidence adapter.

Validates, normalizes, and translates a gene-level exome association file
into the --exomes-in argument for priors.py.

Engine contract (traced from priors.py calculate_huge_scores_exomes, line 4930):
  - Header required
  - Gene ID column: auto-detected or --exomes-gene-col
  - Requires at least 2 of: p-value, beta, SE (the third is derived internally)
  - Auto-detects column names via _determine_columns() if not specified
  - Optional: --exomes-n or --exomes-n-col for sample size (SE = 1/sqrt(N))
  - Duplicate genes: first occurrence kept, subsequent skipped with warning
  - Processing: p/beta/SE → Z → log-BF (HuGE score) internally
  - Whitespace-delimited
  - "NA" values silently skipped
  - p-values: clamp < 1e-250, reject ≤ 0 or > 1
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


class ExomeAdapter(BaseEvidenceAdapter):
    """Adapter for gene-level exome/rare-variant association inputs."""

    ANALYSIS_TYPE = "exome"
    ENGINE_FLAG = "--exomes-in"
    EVIDENCE_DESCRIPTION = (
        "Gene-level exome association statistics (gene, p-value, beta, SE). "
        "Engine computes HuGE scores internally. Requires at least 2 of: "
        "p-value, beta, SE columns."
    )

    def validate_file(self, input_path, params=None):
        """Validate the exome evidence file.

        The engine auto-detects columns, so we don't require specific column
        names — but we verify the file is a valid tabular file with a header
        and at least one data row.
        """
        params = params or {}

        result = validate_tabular_file(input_path, min_data_rows=1)
        if result.status == ValidationStatus.FAIL:
            return result

        gene_col = params.get("exomes_gene_col")
        p_col = params.get("exomes_p_col")
        beta_col = params.get("exomes_beta_col")
        se_col = params.get("exomes_se_col")

        with open(input_path, "r") as f:
            header = f.readline().strip().split()

        issues = []
        if gene_col and gene_col not in header:
            issues.append(
                f"Specified gene column '{gene_col}' not found in header: "
                f"{', '.join(header)}"
            )
        if p_col and p_col not in header:
            issues.append(
                f"Specified p-value column '{p_col}' not found in header: "
                f"{', '.join(header)}"
            )
        if beta_col and beta_col not in header:
            issues.append(
                f"Specified beta column '{beta_col}' not found in header: "
                f"{', '.join(header)}"
            )
        if se_col and se_col not in header:
            issues.append(
                f"Specified SE column '{se_col}' not found in header: "
                f"{', '.join(header)}"
            )

        if issues:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=issues,
            )

        return FileValidationResult(status=ValidationStatus.PASS)

    def normalize(self, input_path, output_dir, params=None):
        """Normalize the exome evidence file.

        Preserves original, writes a clean copy with blank lines removed.
        The engine handles all column detection and statistical computation
        internally — the wrapper does NOT transform values.
        """
        params = params or {}

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, "exome_original.tsv")
        shutil.copy2(input_path, original_dest)

        warnings = []
        transformations = [
            "No wrapper-level transformation. Engine computes HuGE scores: "
            "p/beta/SE → Z → log-BF internally."
        ]

        with open(input_path, "r") as f:
            all_lines = f.readlines()

        if not all_lines:
            return EvidenceNormalizationResult(
                original_path=original_dest,
                normalized_path=os.path.join(norm_dir, "exome.tsv"),
                original_rows=0, normalized_rows=0, removed_rows=0,
                transformations=transformations, warnings=["Empty file"],
            )

        header = all_lines[0]
        data_lines = all_lines[1:]
        original_rows = len(data_lines)
        blanks = 0
        valid_lines = []

        for line in data_lines:
            if not line.strip():
                blanks += 1
                continue
            valid_lines.append(line)

        normalized_path = os.path.join(norm_dir, "exome.tsv")
        with open(normalized_path, "w") as f:
            f.write(header)
            for line in valid_lines:
                f.write(line)

        return EvidenceNormalizationResult(
            original_path=original_dest,
            normalized_path=normalized_path,
            original_rows=original_rows,
            normalized_rows=len(valid_lines),
            removed_rows=blanks,
            blanks_removed=blanks,
            transformations=transformations,
            warnings=warnings,
        )

    def build_engine_args(self, normalized_path, params=None):
        """Return priors.py CLI arguments for exome input.

        Passes the file via --exomes-in. Optional column overrides are
        passed through if specified.
        """
        params = params or {}
        args = ["--exomes-in", normalized_path]

        if params.get("exomes_gene_col"):
            args.extend(["--exomes-gene-col", params["exomes_gene_col"]])
        if params.get("exomes_p_col"):
            args.extend(["--exomes-p-col", params["exomes_p_col"]])
        if params.get("exomes_beta_col"):
            args.extend(["--exomes-beta-col", params["exomes_beta_col"]])
        if params.get("exomes_se_col"):
            args.extend(["--exomes-se-col", params["exomes_se_col"]])
        if params.get("exomes_n_col"):
            args.extend(["--exomes-n-col", params["exomes_n_col"]])
        if params.get("exomes_n") is not None:
            args.extend(["--exomes-n", str(params["exomes_n"])])

        return args
