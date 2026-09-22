"""Input gene quality control for the PIGEAN wrapper.

Validates normalized genes against the gene map, location files,
and annotation gene sets to produce a per-gene QC report.
"""

import csv
import os

from pigean.adapters import GeneQCResult, BackgroundQCResult, ValidationStatus


def load_gene_map(gene_map_path):
    """Load known gene symbols from the gene map file.

    Format: source_id<tab>gene_symbol (2 columns).
    Returns the set of unique gene symbols (column 2).
    """
    symbols = set()
    with open(gene_map_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                symbols.add(parts[1])
    return symbols


def load_gene_locations(loc_file_path):
    """Load gene symbols that have genomic coordinates.

    Format: gene_id<tab>chrom<tab>start<tab>end<tab>strand<tab>gene_symbol
    (6 columns, gene symbol in column 6).
    Returns the set of gene symbols with coordinates.
    """
    symbols = set()
    with open(loc_file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                symbols.add(parts[5])
    return symbols


def count_gene_set_memberships(genes, gene_set_files):
    """Count how many gene sets contain each gene.

    Each gene-set file has lines: gene_set_name<tab>gene1<tab>gene2<tab>...
    Returns dict mapping gene → count of gene sets containing it.
    """
    counts = {gene: 0 for gene in genes}
    gene_set = set(genes)

    for gs_file in gene_set_files:
        with open(gs_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                # First column is gene set name; remaining are genes
                members = set(parts[1:])
                for gene in gene_set & members:
                    counts[gene] += 1

    return counts


def validate_input_genes(normalized_path, gene_map_path, gene_loc_path,
                         gene_set_files):
    """Run full QC on normalized input genes.

    Checks each gene for:
    1. Recognition in the gene map
    2. Coordinates in the gene location file
    3. Membership in annotation gene sets

    Returns GeneQCResult with per-gene details and aggregate counts.
    """
    # Load normalized genes
    with open(normalized_path, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    if not genes:
        return GeneQCResult(
            status=ValidationStatus.FAIL,
            issues=["No genes in normalized input"],
            input_count=0,
        )

    # Load references
    known_symbols = load_gene_map(gene_map_path)
    gene_locs = load_gene_locations(gene_loc_path)
    memberships = count_gene_set_memberships(genes, gene_set_files)

    per_gene = []
    recognized_count = 0
    with_coordinates_count = 0
    with_annotations_count = 0
    issues = []

    for gene in genes:
        recognized = gene in known_symbols
        has_coords = gene in gene_locs
        n_memberships = memberships.get(gene, 0)

        if recognized:
            recognized_count += 1
        if has_coords:
            with_coordinates_count += 1
        if n_memberships > 0:
            with_annotations_count += 1

        # Determine per-gene status and notes
        notes = []
        if not recognized:
            status = "UNRESOLVED"
            notes.append("not found in gene map")
        elif not has_coords:
            status = "WARNING"
            notes.append("no genomic coordinates")
        elif n_memberships == 0:
            status = "WARNING"
            notes.append("not in any annotation gene set")
        else:
            status = "OK"

        per_gene.append({
            "gene": gene,
            "recognized": recognized,
            "has_coordinates": has_coords,
            "annotation_memberships": n_memberships,
            "status": status,
            "notes": "; ".join(notes) if notes else "",
        })

    unresolved_count = len(genes) - recognized_count

    # Determine overall status
    if recognized_count == 0:
        overall_status = ValidationStatus.FAIL
        issues.append("No recognized genes in input")
    elif unresolved_count > 0:
        overall_status = ValidationStatus.PASS_WITH_WARNINGS
        unresolved_names = [g for g in genes if g not in known_symbols]
        issues.append(
            f"{unresolved_count} gene(s) not recognized: "
            + ", ".join(unresolved_names)
        )
    else:
        overall_status = ValidationStatus.PASS

    return GeneQCResult(
        status=overall_status,
        issues=issues,
        input_count=len(genes),
        recognized_count=recognized_count,
        unresolved_count=unresolved_count,
        with_coordinates_count=with_coordinates_count,
        with_annotations_count=with_annotations_count,
        per_gene=per_gene,
    )


def validate_background(background_normalized_path, input_normalized_path,
                         gene_map_path):
    """Check background gene consistency with input genes.

    Validates:
    1. How many background genes are recognized in the gene map.
    2. Whether all input (positive-control) genes are present in the background.

    Parameters
    ----------
    background_normalized_path : str
        Path to the normalized background gene list.
    input_normalized_path : str
        Path to the normalized input (positive-control) gene list.
    gene_map_path : str
        Path to the gene map file.

    Returns
    -------
    BackgroundQCResult
    """
    # Load background genes
    with open(background_normalized_path, "r") as f:
        bg_genes = [line.strip() for line in f if line.strip()]

    # Load input genes
    with open(input_normalized_path, "r") as f:
        input_genes = [line.strip() for line in f if line.strip()]

    # Load gene map
    known_symbols = load_gene_map(gene_map_path)

    bg_set = set(bg_genes)
    total_count = len(bg_genes)
    mapped_count = sum(1 for g in bg_genes if g in known_symbols)
    unmapped_count = total_count - mapped_count

    # Check which input genes are missing from the background
    missing = [g for g in input_genes if g not in bg_set]

    issues = []
    if mapped_count == 0:
        status = ValidationStatus.FAIL
        issues.append("No recognized genes in background file")
    elif missing:
        status = ValidationStatus.PASS_WITH_WARNINGS
        issues.append(
            f"{len(missing)} input gene(s) not found in background: "
            + ", ".join(missing)
        )
    else:
        status = ValidationStatus.PASS

    return BackgroundQCResult(
        status=status,
        issues=issues,
        total_count=total_count,
        mapped_count=mapped_count,
        unmapped_count=unmapped_count,
        input_genes_missing_from_background=missing,
    )


def validate_evidence_genes(normalized_path, gene_map_path, config):
    """Run gene QC on a tabular evidence file (BF, Z-score, percentile, exome).

    Checks each gene for recognition in the gene map.
    Returns GeneQCResult with per-gene details.
    """
    gene_col = config.get("gene_column", "Gene")
    analysis = config.get("analysis", "")

    if analysis == "exome":
        gene_col = config.get("exomes_gene_col") or gene_col

    known_symbols = load_gene_map(gene_map_path)

    with open(normalized_path, "r") as f:
        header = f.readline().strip().split()
        lines = f.readlines()

    if gene_col not in header:
        for col in header:
            if col.lower() in ("gene", "gene_id", "gene_symbol", "geneid"):
                gene_col = col
                break

    if gene_col not in header:
        return GeneQCResult(
            status=ValidationStatus.PASS_WITH_WARNINGS,
            issues=[f"Could not identify gene column in evidence file (tried '{gene_col}')"],
            input_count=len(lines),
        )

    gene_idx = header.index(gene_col)
    seen = set()
    per_gene = []
    recognized_count = 0
    duplicate_count = 0

    for line in lines:
        cols = line.strip().split()
        if not cols or gene_idx >= len(cols):
            continue
        gene = cols[gene_idx]
        if gene in seen:
            duplicate_count += 1
            continue
        seen.add(gene)
        recognized = gene in known_symbols
        if recognized:
            recognized_count += 1
        per_gene.append({
            "gene": gene,
            "recognized": recognized,
            "has_coordinates": False,
            "annotation_memberships": 0,
            "status": "OK" if recognized else "UNRESOLVED",
            "notes": "" if recognized else "not found in gene map",
        })

    input_count = len(per_gene)
    unresolved_count = input_count - recognized_count
    issues = []

    if input_count == 0:
        return GeneQCResult(
            status=ValidationStatus.FAIL,
            issues=["No genes found in evidence file"],
            input_count=0,
        )

    if recognized_count == 0:
        return GeneQCResult(
            status=ValidationStatus.FAIL,
            issues=["No recognized genes in evidence file"],
            input_count=input_count,
            per_gene=per_gene,
        )

    if unresolved_count > 0:
        unresolved_names = [g["gene"] for g in per_gene if g["status"] == "UNRESOLVED"]
        issues.append(
            f"{unresolved_count} gene(s) not recognized: "
            + ", ".join(unresolved_names[:20])
            + ("..." if unresolved_count > 20 else "")
        )

    return GeneQCResult(
        status=(ValidationStatus.PASS_WITH_WARNINGS if issues
                else ValidationStatus.PASS),
        issues=issues,
        input_count=input_count,
        recognized_count=recognized_count,
        unresolved_count=unresolved_count,
        per_gene=per_gene,
    )


def write_evidence_qc_tsv(qc_result, output_dir):
    """Write input_evidence_qc.tsv to the output directory."""
    path = os.path.join(output_dir, "input_evidence_qc.tsv")
    fieldnames = [
        "gene", "recognized", "status", "notes",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in qc_result.per_gene:
            writer.writerow({k: row.get(k, "") for k in fieldnames})
    return path


def write_gene_qc_tsv(qc_result, output_dir):
    """Write input_gene_qc.tsv to the output directory.

    Columns: gene, recognized, has_coordinates, annotation_memberships,
             status, notes
    """
    path = os.path.join(output_dir, "input_gene_qc.tsv")
    fieldnames = [
        "gene", "recognized", "has_coordinates",
        "annotation_memberships", "status", "notes",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in qc_result.per_gene:
            writer.writerow(row)
    return path


def compute_validation_status(file_result, qc_result,
                               background_file_result=None,
                               background_qc_result=None,
                               custom_gs_prep_result=None):
    """Combine all validation stages into final ValidationStatus.

    Parameters
    ----------
    file_result : FileValidationResult
        Primary input file validation.
    qc_result : GeneQCResult
        Gene QC result.
    background_file_result : FileValidationResult or None
        Background file validation (if background was provided).
    background_qc_result : BackgroundQCResult or None
        Background consistency check result.
    custom_gs_prep_result : GeneSetPreparationResult or None
        Custom gene-set preparation result.

    Returns
    -------
    ValidationStatus
    """
    # Check for FAIL in any stage
    if file_result.status == ValidationStatus.FAIL:
        return ValidationStatus.FAIL

    if qc_result.status == ValidationStatus.FAIL:
        return ValidationStatus.FAIL

    if (background_file_result is not None
            and background_file_result.status == ValidationStatus.FAIL):
        return ValidationStatus.FAIL

    if (background_qc_result is not None
            and background_qc_result.status == ValidationStatus.FAIL):
        return ValidationStatus.FAIL

    if (custom_gs_prep_result is not None
            and custom_gs_prep_result.status == ValidationStatus.FAIL):
        return ValidationStatus.FAIL

    # Check for warnings in any stage
    if (file_result.status == ValidationStatus.PASS_WITH_WARNINGS
            or qc_result.status == ValidationStatus.PASS_WITH_WARNINGS):
        return ValidationStatus.PASS_WITH_WARNINGS

    if (background_qc_result is not None
            and background_qc_result.status == ValidationStatus.PASS_WITH_WARNINGS):
        return ValidationStatus.PASS_WITH_WARNINGS

    if (custom_gs_prep_result is not None
            and custom_gs_prep_result.status == ValidationStatus.PASS_WITH_WARNINGS):
        return ValidationStatus.PASS_WITH_WARNINGS

    return ValidationStatus.PASS
