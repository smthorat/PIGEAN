"""Positive-control (gene list) evidence adapter.

Validates, normalizes, and translates a simple gene-per-line input file
into the --positive-controls-in argument for priors.py.
"""

import os
import shutil

from pigean.adapters import (
    FileValidationResult,
    NormalizationResult,
    ValidationStatus,
)
from pigean.adapters.base import BaseEvidenceAdapter, EvidenceNormalizationResult


class PositiveControlsAdapter(BaseEvidenceAdapter):
    """Adapter for positive-control gene list inputs."""

    ANALYSIS_TYPE = "positive-controls"
    ENGINE_FLAG = "--positive-controls-in"
    EVIDENCE_DESCRIPTION = "Gene list (one gene symbol per line)"

    def validate_file(self, input_path, params=None):
        """Validate the input file at the file level.

        Checks:
        - path exists
        - path is a regular file
        - file is readable
        - file is non-empty
        - file appears to be text (no null bytes in first 8 KB)
        - file has one-gene-per-line structure (short lines, no tabs)

        Returns FileValidationResult with PASS or FAIL.
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

        file_size = os.path.getsize(input_path)
        if file_size == 0:
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

        non_blank_lines = [ln for ln in lines if ln.strip()]
        if not non_blank_lines:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=["Input file contains only blank lines"],
            )

        for i, line in enumerate(non_blank_lines, 1):
            stripped = line.strip()
            if "\t" in stripped:
                issues.append(
                    f"Line {i} contains tab characters (expected one gene per line)"
                )
            if len(stripped) > 100:
                issues.append(
                    f"Line {i} is unusually long ({len(stripped)} chars)"
                )

        if issues:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=issues,
            )

        return FileValidationResult(status=ValidationStatus.PASS)

    def normalize(self, input_path, output_dir, params=None,
                  input_label="positive_controls",
                  normalized_label="positive_controls"):
        """Normalize the gene list input.

        1. Copies original to output_dir/input/<input_label>.txt
        2. Strips whitespace, removes blank lines, deduplicates (preserving
           first-occurrence order)
        3. Writes to output_dir/normalized/<normalized_label>.txt

        Parameters
        ----------
        input_label : str
            Filename stem for the preserved original. Default "positive_controls".
            Use "background_original" for background files.
        normalized_label : str
            Filename stem for the normalized output. Default "positive_controls".
            Use "background" for background files.

        Returns NormalizationResult.
        """
        warnings = []

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, f"{input_label}.txt")
        shutil.copy2(input_path, original_dest)

        with open(input_path, "r") as f:
            raw_lines = f.readlines()

        original_count = len(raw_lines)
        blanks_removed = 0
        duplicates_removed = 0
        seen = set()
        normalized_genes = []

        for line in raw_lines:
            gene = line.strip()
            if not gene:
                blanks_removed += 1
                continue
            if gene in seen:
                duplicates_removed += 1
                warnings.append(f"Duplicate gene removed: {gene}")
                continue
            seen.add(gene)
            normalized_genes.append(gene)

        if blanks_removed > 0:
            warnings.append(f"Removed {blanks_removed} blank line(s)")
        if duplicates_removed > 0:
            warnings.append(
                f"Removed {duplicates_removed} duplicate gene(s)"
            )

        normalized_path = os.path.join(norm_dir, f"{normalized_label}.txt")
        with open(normalized_path, "w") as f:
            for gene in normalized_genes:
                f.write(gene + "\n")

        return NormalizationResult(
            normalized_path=normalized_path,
            original_count=original_count,
            normalized_count=len(normalized_genes),
            duplicates_removed=duplicates_removed,
            blanks_removed=blanks_removed,
            warnings=warnings,
        )

    def build_engine_args(self, normalized_path, params=None):
        """Return the priors.py CLI arguments for positive-control input.

        Returns ["--positive-controls-in", <normalized_path>].
        """
        return ["--positive-controls-in", normalized_path]
