"""Run provenance manifest for the PIGEAN wrapper.

Every attempted run — including validation failures where the engine
was never executed — must produce a run_manifest.json.
"""

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone

import pigean


def compute_file_sha256(filepath):
    """Compute the SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def get_dependency_versions():
    """Return a dict of Python/numpy/scipy versions."""
    versions = {
        "python": platform.python_version(),
    }
    try:
        import numpy
        versions["numpy"] = numpy.__version__
    except ImportError:
        versions["numpy"] = "NOT_INSTALLED"
    try:
        import scipy
        versions["scipy"] = scipy.__version__
    except ImportError:
        versions["scipy"] = "NOT_INSTALLED"
    return versions


def generate_manifest(config, file_result, norm_result, qc_result,
                      validation_status, validation_issues,
                      execution_status, exit_code,
                      base_dir, reference_paths, gene_set_paths,
                      gene_map_path, output_dir,
                      start_time, end_time,
                      background_path=None,
                      background_norm_result=None,
                      background_qc_result=None,
                      custom_gs_prep_result=None,
                      evidence_norm_result=None,
                      parsed_outputs=None,
                      convergence_result=None,
                      interpretation_result=None,
                      advanced_mode=None):
    """Build the full provenance manifest dict.

    Parameters
    ----------
    config : dict
        Resolved configuration.
    file_result : FileValidationResult or None
        File-level validation result.
    norm_result : NormalizationResult or None
        Normalization result.
    qc_result : GeneQCResult or None
        Gene QC result.
    validation_status : str
        Overall validation status (PASS / PASS_WITH_WARNINGS / FAIL).
    validation_issues : list of str
        All validation issues.
    execution_status : str
        NOT_RUN / COMPLETED / FAILED.
    exit_code : int or None
        Engine process exit code (None if not run).
    base_dir : str
        Application base directory.
    reference_paths : dict or None
        Resolved reference paths.
    gene_set_paths : list or None
        Resolved gene-set paths.
    gene_map_path : str or None
        Resolved gene map path.
    output_dir : str
        Output directory.
    start_time : datetime
        Run start time.
    end_time : datetime
        Run end time.
    background_path : str or None
        Original background gene list file path.
    background_norm_result : NormalizationResult or None
        Background normalization result.
    background_qc_result : BackgroundQCResult or None
        Background consistency check result.
    custom_gs_prep_result : GeneSetPreparationResult or None
        Custom gene-set preparation result.

    Returns
    -------
    dict
        Complete manifest.
    """
    priors_path = os.path.join(base_dir, "engine", "priors.py")

    manifest = {
        "wrapper_version": pigean.__version__,
        "engine": {
            "name": "priors.py",
            "sha256": (compute_file_sha256(priors_path)
                       if os.path.isfile(priors_path) else "UNAVAILABLE"),
        },
        "dependencies": get_dependency_versions(),
        "analysis": {
            "analysis_type": config.get("analysis", "positive-controls"),
            "model_mode": config.get("mode", "standard"),
            "preset": config.get("preset", "standard"),
            "genome_build": config.get("genome_build", "hg19"),
            "gene_sets": config.get("gene_sets", "default"),
        },
        "settings": {
            "max_num_gene_sets": config.get("max_num_gene_sets"),
            "gene_filter_value": config.get("gene_filter_value"),
            "gene_set_filter_value": config.get("gene_set_filter_value"),
            "debug_level": config.get("debug_level"),
        },
        "input": {},
        "references": {},
        "gene_qc": {},
        "validation": {
            "status": (validation_status.value
                       if hasattr(validation_status, "value")
                       else str(validation_status)),
            "issues": validation_issues,
        },
        "execution": {
            "status": execution_status,
            "exit_code": exit_code,
        },
        "timing": {
            "start": start_time.isoformat() if start_time else None,
            "end": end_time.isoformat() if end_time else None,
            "duration_seconds": (
                round((end_time - start_time).total_seconds(), 2)
                if start_time and end_time else None
            ),
        },
        "log_file": "pigean_run.log",
        "stability_status": (
            convergence_result["stability_status"].value
            if convergence_result is not None
            and hasattr(convergence_result.get("stability_status"), "value")
            else (str(convergence_result["stability_status"])
                  if convergence_result is not None
                  else "NOT_ASSESSED")
        ),
        "formal_convergence_status": (
            convergence_result["formal_convergence_status"].value
            if convergence_result is not None
            and hasattr(convergence_result.get("formal_convergence_status"),
                        "value")
            else "NOT_FORMALLY_ASSESSED"
        ),
    }

    # Phase 6: high-level model mode provenance. This remains separate from
    # analysis/evidence type and records only verified mode-level settings.
    if advanced_mode is None:
        advanced_mode = {
            "name": config.get("mode", "standard"),
            "user_supplied_parameters": {},
            "resolved_parameters": {},
            "engine_arguments": [],
            "compatibility_status": "UNKNOWN",
            "stability_applicability": "UNKNOWN",
            "interpretation_compatibility": "UNKNOWN",
        }
    manifest["advanced_mode"] = advanced_mode

    # Input checksums
    if norm_result:
        original_path = os.path.join(output_dir, "input",
                                     "positive_controls.txt")
        if os.path.isfile(original_path):
            manifest["input"]["original_checksum"] = (
                compute_file_sha256(original_path)
            )
        if os.path.isfile(norm_result.normalized_path):
            manifest["input"]["normalized_checksum"] = (
                compute_file_sha256(norm_result.normalized_path)
            )
        manifest["input"]["original_count"] = norm_result.original_count
        manifest["input"]["normalized_count"] = norm_result.normalized_count
        manifest["input"]["duplicates_removed"] = norm_result.duplicates_removed
        manifest["input"]["blanks_removed"] = norm_result.blanks_removed

    # Reference checksums
    if gene_map_path and os.path.isfile(gene_map_path):
        manifest["references"]["gene_map_checksum"] = (
            compute_file_sha256(gene_map_path)
        )
    if gene_set_paths:
        manifest["references"]["gene_set_checksums"] = {
            os.path.basename(p): compute_file_sha256(p)
            for p in gene_set_paths if os.path.isfile(p)
        }
    if reference_paths:
        for key, path in reference_paths.items():
            if os.path.isfile(path):
                manifest["references"][f"{key}_checksum"] = (
                    compute_file_sha256(path)
                )

    # Gene QC summary — handle both GeneQCResult and EvidenceQCResult
    if qc_result:
        if hasattr(qc_result, "input_count"):
            # GeneQCResult (positive-controls)
            manifest["gene_qc"] = {
                "input_count": qc_result.input_count,
                "recognized_count": qc_result.recognized_count,
                "unresolved_count": qc_result.unresolved_count,
                "with_coordinates_count": qc_result.with_coordinates_count,
                "with_annotations_count": qc_result.with_annotations_count,
            }
        elif hasattr(qc_result, "total_genes"):
            # EvidenceQCResult (BF, Z-score, percentile, exome, GWAS)
            manifest["gene_qc"] = {
                "evidence_type": getattr(qc_result, "evidence_type", "unknown"),
                "total_genes": qc_result.total_genes,
                "recognized_genes": qc_result.recognized_genes,
                "unresolved_genes": qc_result.unresolved_genes,
                "duplicate_genes": qc_result.duplicate_genes,
                "missing_values": qc_result.missing_values,
                "invalid_values": qc_result.invalid_values,
            }

    # Custom gene-set provenance (Phase 2)
    if config.get("custom_gene_set_format"):
        custom_gs = {
            "format": config.get("custom_gene_set_format"),
            "action": config.get("custom_gene_set_action"),
            "file_count": (len(custom_gs_prep_result.original_paths)
                           if custom_gs_prep_result else 0),
            "files": {},
        }
        if custom_gs_prep_result and custom_gs_prep_result.original_paths:
            for p in custom_gs_prep_result.original_paths:
                if os.path.isfile(p):
                    custom_gs["files"][os.path.basename(p)] = (
                        compute_file_sha256(p)
                    )
        manifest["custom_gene_sets"] = custom_gs

    # Background provenance (Phase 2)
    if background_path is not None:
        bg_section = {
            "provided": True,
            "original_checksum": (
                compute_file_sha256(background_path)
                if os.path.isfile(background_path) else None
            ),
            "normalized_checksum": (
                compute_file_sha256(background_norm_result.normalized_path)
                if (background_norm_result
                    and os.path.isfile(background_norm_result.normalized_path))
                else None
            ),
        }
        if background_qc_result:
            bg_section["qc"] = {
                "total_count": background_qc_result.total_count,
                "mapped_count": background_qc_result.mapped_count,
                "unmapped_count": background_qc_result.unmapped_count,
                "input_genes_missing": (
                    background_qc_result.input_genes_missing_from_background
                ),
            }
        manifest["background"] = bg_section

    # Evidence provenance (Phase 3)
    if evidence_norm_result is not None:
        evidence_section = {
            "type": config.get("analysis", "unknown"),
            "original_file": evidence_norm_result.original_path,
            "original_checksum": (
                compute_file_sha256(evidence_norm_result.original_path)
                if os.path.isfile(evidence_norm_result.original_path)
                else None
            ),
            "normalized_file": evidence_norm_result.normalized_path,
            "normalized_checksum": (
                compute_file_sha256(evidence_norm_result.normalized_path)
                if os.path.isfile(evidence_norm_result.normalized_path)
                else None
            ),
            "original_rows": evidence_norm_result.original_rows,
            "normalized_rows": evidence_norm_result.normalized_rows,
            "removed_rows": evidence_norm_result.removed_rows,
            "duplicates_removed": evidence_norm_result.duplicates_removed,
            "transformations": evidence_norm_result.transformations,
            "warnings": evidence_norm_result.warnings,
        }
        if config.get("gene_column"):
            evidence_section["gene_column"] = config["gene_column"]
        if config.get("score_column"):
            evidence_section["score_column"] = config["score_column"]
        manifest["evidence"] = evidence_section

    # GWAS column override provenance (Phase 4)
    if config.get("analysis") == "gwas":
        gwas_section = {}
        for key in ["gwas_chrom_col", "gwas_pos_col", "gwas_p_col",
                     "gwas_beta_col", "gwas_se_col", "gwas_n_col",
                     "gwas_freq_col", "gwas_locus_col",
                     "gwas_filter_col", "gwas_filter_value"]:
            val = config.get(key)
            if val:
                gwas_section[key] = val
        gwas_n = config.get("gwas_n")
        if gwas_n is not None:
            gwas_section["gwas_n"] = gwas_n
        if gwas_section:
            manifest["gwas_column_overrides"] = gwas_section

    # Phase 5: Stability assessment (post-audit terminology)
    if convergence_result is not None:
        manifest["engine_stability"] = {
            "status": (convergence_result["stability_status"].value
                       if hasattr(convergence_result["stability_status"], "value")
                       else str(convergence_result["stability_status"])),
            "method": convergence_result.get("stability_method",
                                             "PIGEAN_MAX_FRACTIONAL_SEM"),
            "assessment_source": convergence_result["assessment_source"],
            "max_fractional_sem": (
                convergence_result.get("details", {}).get("final_max_sem_ratio")
            ),
            "threshold": convergence_result.get("details", {}).get(
                "max_frac_sem_threshold", 0.01),
            "details": convergence_result.get("details", {}),
            "evidence": convergence_result.get("evidence", []),
            "warnings": convergence_result.get("warnings", []),
        }
        manifest["formal_mcmc"] = {
            "status": convergence_result["formal_convergence_status"].value,
            "reason": convergence_result.get("formal_convergence_reason"),
        }
        manifest["iterations"] = {
            "observed": convergence_result.get("details", {}).get(
                "num_gibbs_iter"),
            "maximum": convergence_result.get("details", {}).get(
                "max_num_iter"),
            "status": convergence_result["iteration_status"].value,
        }
        # Update top-level status keys
        manifest["stability_status"] = (
            manifest["engine_stability"]["status"])
        manifest["formal_convergence_status"] = (
            manifest["formal_mcmc"]["status"])

    # Phase 5: Output statistics
    if parsed_outputs is not None:
        output_stats = {}
        if parsed_outputs.get("gs") is not None:
            gs = parsed_outputs["gs"]
            output_stats["gs_row_count"] = gs.get("row_count", 0)
            output_stats["gs_columns"] = gs.get("columns", [])
            output_stats["has_positive_control_column"] = gs.get(
                "has_positive_control", False)
        if parsed_outputs.get("gss") is not None:
            gss = parsed_outputs["gss"]
            output_stats["gss_row_count"] = gss.get("row_count", 0)
        if parsed_outputs.get("ggss") is not None:
            ggss = parsed_outputs["ggss"]
            output_stats["ggss_row_count"] = ggss.get("row_count", 0)
        if parsed_outputs.get("p") is not None:
            p = parsed_outputs["p"]
            output_stats["p_param_count"] = len(p.get("params", {}))
        manifest["output_statistics"] = output_stats

    # Phase 5: Interpretation summary
    if interpretation_result is not None:
        gs_summary = interpretation_result.get("gene_summary")
        interp_summary = {
            "analysis_type": interpretation_result.get("analysis_type"),
            "convergence_note": interpretation_result.get("convergence_note"),
        }
        if gs_summary is not None:
            interp_summary["top_input_genes_count"] = len(
                gs_summary.get("input_genes", []))
            interp_summary["top_candidate_genes_count"] = len(
                gs_summary.get("candidate_genes", []))
            interp_summary["total_genes"] = gs_summary.get("total_genes", 0)
        gss_summary = interpretation_result.get("gene_set_summary")
        if gss_summary is not None:
            interp_summary["top_gene_sets_count"] = len(
                gss_summary.get("top_gene_sets", []))
            interp_summary["total_gene_sets"] = gss_summary.get(
                "total_gene_sets", 0)
        manifest["interpretation_summary"] = interp_summary

    return manifest


def write_manifest(manifest, output_dir):
    """Write run_manifest.json to the output directory."""
    path = os.path.join(output_dir, "run_manifest.json")
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
    return path
