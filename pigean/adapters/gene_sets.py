"""Custom gene-set validation, conversion, and preparation.

Handles GMT→engine format conversion and enforces the legal
combinations table for custom gene-set configuration.
"""

import os
import shutil

from pigean.adapters import (
    FileValidationResult,
    GeneSetPreparationResult,
    ValidationStatus,
)


def validate_gene_set_file(file_path):
    """Validate a single gene-set file exists, is readable, non-empty, and text.

    Unlike PositiveControlsAdapter.validate_file, does NOT check for
    one-gene-per-line structure (gene-set files have tabs).

    Parameters
    ----------
    file_path : str
        Absolute path to a gene-set file.

    Returns
    -------
    FileValidationResult
        PASS if the file is valid, FAIL with issues otherwise.
    """
    if not os.path.exists(file_path):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Gene-set file does not exist: {file_path}"],
        )

    if not os.path.isfile(file_path):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Gene-set path is not a regular file: {file_path}"],
        )

    if not os.access(file_path, os.R_OK):
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Gene-set file is not readable: {file_path}"],
        )

    if os.path.getsize(file_path) == 0:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Gene-set file is empty: {file_path}"],
        )

    # Check for binary content (null bytes in first 8 KB)
    with open(file_path, "rb") as f:
        chunk = f.read(8192)
    if b"\x00" in chunk:
        return FileValidationResult(
            status=ValidationStatus.FAIL,
            issues=[f"Gene-set file appears to be binary: {file_path}"],
        )

    return FileValidationResult(status=ValidationStatus.PASS)


def convert_gmt_to_engine_format(gmt_path, output_path):
    """Convert GMT format to engine format by stripping the description column.

    GMT format:  gene_set_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...
    Engine format: gene_set_name<TAB>gene1<TAB>gene2<TAB>...

    Parameters
    ----------
    gmt_path : str
        Path to input GMT file.
    output_path : str
        Path to write the converted engine-format file.

    Returns
    -------
    dict
        {"gene_sets_converted": int, "warnings": list}

    Raises
    ------
    ValueError
        If any non-blank line has fewer than 3 tab-separated columns.
    """
    gene_sets_converted = 0
    warnings = []

    with open(gmt_path, "r") as fin, open(output_path, "w") as fout:
        for line_num, line in enumerate(fin, 1):
            stripped = line.strip()
            if not stripped:
                continue

            fields = stripped.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"GMT line {line_num} in {os.path.basename(gmt_path)} has "
                    f"{len(fields)} field(s) — expected at least 3 "
                    f"(name, description, gene1). Line: {stripped[:80]}"
                )

            # Strip description (column index 1)
            name = fields[0]
            genes = fields[2:]
            fout.write(name + "\t" + "\t".join(genes) + "\n")
            gene_sets_converted += 1

    return {"gene_sets_converted": gene_sets_converted, "warnings": warnings}


def prepare_custom_gene_sets(custom_files, custom_format, custom_action,
                              gene_set_profile, output_dir):
    """Validate gene-set configuration and prepare normalized paths.

    Called unconditionally (even with no custom files) so the full
    legal-combinations table is enforced in one place.

    Legal combinations:
        profile=default/mouse-only/msigdb-only + no files + no action → PASS (no-op)
        profile=custom + no files                                     → FAIL
        files provided + no format                                    → FAIL
        files provided + no action                                    → FAIL
        action provided + no files                                    → FAIL
        profile=custom + files + format + action                      → validate → convert → PASS
        profile!=custom + files + format + action                     → validate → convert → PASS

    Parameters
    ----------
    custom_files : list of str
        Paths to custom gene-set files (empty list if none).
    custom_format : str or None
        "engine" or "gmt". Required when custom_files is non-empty.
    custom_action : str or None
        "replace" or "supplement". Required when custom_files is non-empty.
    gene_set_profile : str
        The resolved gene-set profile name.
    output_dir : str
        Run output directory (for writing normalized files).

    Returns
    -------
    GeneSetPreparationResult
    """
    has_files = bool(custom_files)

    # No custom files and not the "custom" profile: no-op
    if not has_files and gene_set_profile != "custom":
        # Also check for orphaned action without files
        if custom_action is not None:
            return GeneSetPreparationResult(
                status=ValidationStatus.FAIL,
                issues=["--custom-gene-set-action requires --custom-gene-set-files"],
            )
        if custom_format is not None:
            return GeneSetPreparationResult(
                status=ValidationStatus.FAIL,
                issues=["--custom-gene-set-format requires --custom-gene-set-files"],
            )
        return GeneSetPreparationResult(status=ValidationStatus.PASS)

    # Custom profile but no files
    if gene_set_profile == "custom" and not has_files:
        return GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=["--gene-sets custom requires --custom-gene-set-files"],
        )

    # Files provided but missing format
    if has_files and custom_format is None:
        return GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=[
                "--custom-gene-set-format is required when "
                "--custom-gene-set-files is provided"
            ],
        )

    # Files provided but missing action
    if has_files and custom_action is None:
        return GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=[
                "--custom-gene-set-action is required when "
                "--custom-gene-set-files is provided"
            ],
        )

    # Action provided but no files (with custom profile)
    if not has_files and custom_action is not None:
        return GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=["--custom-gene-set-action requires --custom-gene-set-files"],
        )

    # Validate each file
    all_issues = []
    for fpath in custom_files:
        result = validate_gene_set_file(fpath)
        if result.status == ValidationStatus.FAIL:
            all_issues.extend(result.issues)

    if all_issues:
        return GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=all_issues,
        )

    # Create directories
    input_dir = os.path.join(output_dir, "input")
    norm_dir = os.path.join(output_dir, "normalized")
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(norm_dir, exist_ok=True)

    # Handle duplicate basenames
    seen_basenames = {}
    original_paths = []
    normalized_paths = []
    gene_sets_converted = 0
    gene_sets_copied = 0
    issues = []

    for fpath in custom_files:
        basename = os.path.basename(fpath)
        stem, ext = os.path.splitext(basename)

        # Deduplicate basenames
        if basename in seen_basenames:
            seen_basenames[basename] += 1
            basename = f"{stem}_{seen_basenames[basename]}{ext}"
            stem = f"{stem}_{seen_basenames[basename]}"
        else:
            seen_basenames[basename] = 0

        # Copy original
        original_dest = os.path.join(input_dir, basename)
        shutil.copy2(fpath, original_dest)
        original_paths.append(original_dest)

        # Convert or copy to normalized
        if custom_format == "gmt":
            norm_basename = f"{stem}.engine.txt"
            norm_dest = os.path.join(norm_dir, norm_basename)
            try:
                result = convert_gmt_to_engine_format(fpath, norm_dest)
                gene_sets_converted += result["gene_sets_converted"]
                issues.extend(result.get("warnings", []))
            except ValueError as e:
                return GeneSetPreparationResult(
                    status=ValidationStatus.FAIL,
                    issues=[str(e)],
                    original_paths=original_paths,
                )
        else:
            # Engine format — copy as-is
            norm_dest = os.path.join(norm_dir, basename)
            shutil.copy2(fpath, norm_dest)
            gene_sets_copied += 1

        normalized_paths.append(norm_dest)

    return GeneSetPreparationResult(
        status=ValidationStatus.PASS,
        original_paths=original_paths,
        normalized_paths=normalized_paths,
        issues=issues,
        gene_sets_converted=gene_sets_converted,
        gene_sets_copied=gene_sets_copied,
    )
