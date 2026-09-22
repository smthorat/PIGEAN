"""GWAS summary statistics evidence adapter.

Validates, normalizes, and translates a SNP-level GWAS summary statistics file
into the --gwas-in argument and associated flags for priors.py.

Engine contract (traced from priors.py calculate_huge_scores_gwas, line 3232):
  - Header required
  - Delimiter: tab if header contains tab, else any whitespace
  - Column auto-detection via _determine_columns() if columns not specified
  - Requires chrom+pos (or locus) AND at least 2 of {p, beta, SE/N}
  - gene_loc_file is MANDATORY (passed separately via reference_paths)
  - Chromosomes cleaned via _clean_chrom() (strips "chr" prefix)
  - Positions: must be integer (cast via int())
  - P-values: clamped >= 1e-250, must be in (0, 1]
  - SE: if missing, can be derived from N as 1/sqrt(N)
  - "NA" values silently skipped
  - Duplicate SNPs at same position: all kept
  - SNP-to-gene mapping: distance-based via gene_loc_file coordinates
  - NO allele or rsID columns required
  - Genome build is NOT validated by the engine — wrapper must enforce
"""

import os
import shutil

from pigean.adapters import FileValidationResult, ValidationStatus
from pigean.adapters.base import (
    BaseEvidenceAdapter,
    EvidenceNormalizationResult,
    validate_tabular_file,
)


class GwasAdapter(BaseEvidenceAdapter):
    """Adapter for GWAS summary statistics (SNP-level) inputs.

    Unlike gene-level adapters (BF, Z-score, percentile, exome), GWAS input
    is SNP-level and requires coordinate-aware reference files (gene_loc_file,
    tss_loc, exons_loc) to map SNPs to genes. The genome build of the GWAS
    file MUST match the reference pack.
    """

    ANALYSIS_TYPE = "gwas"
    ENGINE_FLAG = "--gwas-in"
    EVIDENCE_DESCRIPTION = (
        "GWAS summary statistics (SNP-level: chrom, pos, p-value, beta, SE). "
        "Engine computes HuGE scores via SNP-to-gene mapping using genomic "
        "coordinates. Requires genome-build-matched reference files."
    )

    # Column overrides that can be specified by the user.
    # Maps wrapper param name → engine CLI flag.
    COLUMN_FLAGS = {
        "gwas_chrom_col": "--gwas-chrom-col",
        "gwas_pos_col": "--gwas-pos-col",
        "gwas_p_col": "--gwas-p-col",
        "gwas_beta_col": "--gwas-beta-col",
        "gwas_se_col": "--gwas-se-col",
        "gwas_n_col": "--gwas-n-col",
        "gwas_freq_col": "--gwas-freq-col",
        "gwas_locus_col": "--gwas-locus-col",
        "gwas_filter_col": "--gwas-filter-col",
        "gwas_filter_value": "--gwas-filter-value",
    }

    # Scalar (non-column) flags.
    SCALAR_FLAGS = {
        "gwas_n": "--gwas-n",
    }

    def validate_file(self, input_path, params=None):
        """Validate the GWAS summary statistics file.

        Checks: exists, readable, non-empty, has header, has ≥1 data row.
        If user specified column overrides, verifies those columns exist
        in the header. The engine's own auto-detection handles unspecified
        columns — the wrapper does NOT duplicate that logic.
        """
        params = params or {}

        result = validate_tabular_file(input_path, min_data_rows=1)
        if result.status == ValidationStatus.FAIL:
            return result

        with open(input_path, "r") as f:
            header_line = f.readline().strip()
        # Match engine's delimiter detection: tab if present, else whitespace
        if "\t" in header_line:
            header = header_line.split("\t")
        else:
            header = header_line.split()
        header = [h for h in header if h]

        issues = []
        for param_key, flag_name in self.COLUMN_FLAGS.items():
            col_name = params.get(param_key)
            if col_name and col_name not in header:
                # Also try: the engine allows 1-based integer indices
                try:
                    idx = int(col_name)
                    if idx < 1 or idx > len(header):
                        issues.append(
                            f"Column index {col_name} ({flag_name}) is out of "
                            f"range. Header has {len(header)} columns."
                        )
                except ValueError:
                    issues.append(
                        f"Specified column '{col_name}' ({flag_name}) not "
                        f"found in header: {', '.join(header[:10])}"
                        f"{'...' if len(header) > 10 else ''}"
                    )

        if issues:
            return FileValidationResult(
                status=ValidationStatus.FAIL,
                issues=issues,
            )

        return FileValidationResult(status=ValidationStatus.PASS)

    def normalize(self, input_path, output_dir, params=None):
        """Normalize the GWAS summary statistics file.

        Preserves original, writes a clean copy with blank lines removed.
        The wrapper does NOT transform GWAS data values — the engine handles
        all statistical computation (p/beta/SE → Z → BF → HuGE scores)
        and SNP-to-gene mapping internally.
        """
        params = params or {}

        input_dir = os.path.join(output_dir, "input")
        norm_dir = os.path.join(output_dir, "normalized")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(norm_dir, exist_ok=True)

        original_dest = os.path.join(input_dir, "gwas_original.tsv")
        shutil.copy2(input_path, original_dest)

        warnings = []
        transformations = [
            "No wrapper-level transformation. Engine computes HuGE scores: "
            "SNP-level p/beta/SE → Z → variant BF → SNP-to-gene mapping → "
            "gene-level log-BF internally."
        ]

        with open(input_path, "r") as f:
            all_lines = f.readlines()

        if not all_lines:
            return EvidenceNormalizationResult(
                original_path=original_dest,
                normalized_path=os.path.join(norm_dir, "gwas.tsv"),
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

        normalized_path = os.path.join(norm_dir, "gwas.tsv")
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
        """Return priors.py CLI arguments for GWAS input.

        Passes the file via --gwas-in. Column overrides are passed through
        if specified. Gene location files (--gene-loc-file, etc.) are NOT
        included here — they come from engine.py via reference_paths.
        """
        params = params or {}
        args = ["--gwas-in", normalized_path]

        # Column overrides
        for param_key, flag_name in self.COLUMN_FLAGS.items():
            value = params.get(param_key)
            if value:
                args.extend([flag_name, str(value)])

        # Scalar flags
        for param_key, flag_name in self.SCALAR_FLAGS.items():
            value = params.get(param_key)
            if value is not None:
                args.extend([flag_name, str(value)])

        return args
