"""Base evidence adapter interface and shared result types.

Every evidence adapter (positive controls, Bayes factors, Z-scores,
percentiles, exomes) conforms to this contract. The wrapper selects
the adapter based on --analysis and delegates validation, normalization,
QC, and engine argument construction to it.
"""

import os
from dataclasses import dataclass, field
from typing import List, Optional

from pigean.adapters import ValidationStatus, FileValidationResult


@dataclass
class EvidenceNormalizationResult:
    """Result of evidence file normalization."""
    original_path: str
    normalized_path: str
    original_rows: int = 0
    normalized_rows: int = 0
    removed_rows: int = 0
    duplicates_removed: int = 0
    blanks_removed: int = 0
    transformations: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class EvidenceQCResult:
    """Result of gene-level QC against reference data."""
    evidence_type: str
    status: ValidationStatus
    issues: List[str] = field(default_factory=list)
    total_genes: int = 0
    recognized_genes: int = 0
    unresolved_genes: int = 0
    duplicate_genes: int = 0
    missing_values: int = 0
    invalid_values: int = 0
    per_gene: List[dict] = field(default_factory=list)


class BaseEvidenceAdapter:
    """Contract that every evidence adapter implements.

    Subclasses must define:
      - ANALYSIS_TYPE: str — the --analysis value that selects this adapter
      - ENGINE_FLAG: str — the primary priors.py CLI flag
      - EVIDENCE_DESCRIPTION: str — human-readable description
    """

    ANALYSIS_TYPE: str = NotImplemented
    ENGINE_FLAG: str = NotImplemented
    EVIDENCE_DESCRIPTION: str = NotImplemented

    def validate_file(self, input_path, params):
        """Validate the input file format.

        Returns FileValidationResult with PASS or FAIL.
        """
        raise NotImplementedError

    def normalize(self, input_path, output_dir, params):
        """Normalize input and preserve both original and normalized copies.

        Returns EvidenceNormalizationResult.
        """
        raise NotImplementedError

    def build_engine_args(self, normalized_path, params):
        """Return priors.py CLI arguments for this evidence type.

        Returns list of str.
        """
        raise NotImplementedError

    def describe_contract(self):
        """Return a dict describing this adapter's input contract."""
        return {
            "analysis_type": self.ANALYSIS_TYPE,
            "engine_flag": self.ENGINE_FLAG,
            "description": self.EVIDENCE_DESCRIPTION,
        }


def validate_tabular_file(input_path, required_columns=None,
                          min_data_rows=1):
    """Shared file-level validation for tabular (header + data) evidence files.

    Checks: exists, readable, non-empty, text, has header, has required
    columns, has at least min_data_rows data rows.

    Returns FileValidationResult.
    """
    issues = []

    if not os.path.exists(input_path):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Input file does not exist: {input_path}"],
        )

    if not os.path.isfile(input_path):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Input path is not a regular file: {input_path}"],
        )

    if not os.access(input_path, os.R_OK):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Input file is not readable: {input_path}"],
        )

    if os.path.getsize(input_path) == 0:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=["Input file is empty"],
        )

    with open(input_path, "rb") as f:
        chunk = f.read(8192)
    if b"\x00" in chunk:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=["Input file appears to be binary (contains null bytes)"],
        )

    with open(input_path, "r") as f:
        lines = f.readlines()

    non_blank = [ln for ln in lines if ln.strip()]
    if not non_blank:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=["Input file contains only blank lines"],
        )

    header_cols = non_blank[0].strip().split()

    if required_columns:
        missing = [c for c in required_columns if c not in header_cols]
        if missing:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=[
                    f"Missing required column(s): {', '.join(missing)}. "
                    f"Header has: {', '.join(header_cols)}. "
                    f"Specify column names with --gene-column / --score-column."
                ],
            )

    data_rows = len(non_blank) - 1
    if data_rows < min_data_rows:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"File has {data_rows} data row(s), need at least {min_data_rows}"],
        )

    return FileValidationResult(status=ValidationStatus.PASS)
