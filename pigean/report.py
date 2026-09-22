"""
pigean/report.py — Text and HTML report generation.

Phase 5: Upgraded with interpretation sections, PIGEAN stability assessment,
top genes/gene-sets, gene-pathway links, metric guide, and caveats.
Only VERIFIED metrics are interpreted. Unverified metrics are not
displayed as interpreted results.

Terminology note (post-audit): The engine's max-fractional-SEM criterion
is reported as "PIGEAN stability" — NOT as formal MCMC convergence.
"""

import os
from datetime import datetime


def generate_report(config, qc_result, manifest, output_dir,
                    evidence_norm_result=None,
                    interpretation_result=None,
                    convergence_result=None):
    """Generate report.txt with QC, execution, interpretation, and stability."""
    lines = []

    lines.append("=" * 60)
    lines.append("  PIGEAN Run Report")
    lines.append("=" * 60)
    lines.append("")

    # Timestamp
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Wrapper version: {manifest.get('wrapper_version', 'unknown')}")
    lines.append("")

    # Run status
    lines.append("RUN STATUS")
    lines.append("-" * 40)
    exec_status = manifest.get("execution", {}).get("status", "UNKNOWN")
    val_status = manifest.get("validation", {}).get("status", "UNKNOWN")

    if convergence_result is not None:
        stab_status = convergence_result["stability_status"].value
        stab_source = convergence_result["assessment_source"]
    else:
        stab_status = "NOT ASSESSED"
        stab_source = None

    lines.append(f"  Validation:  {val_status}")
    lines.append(f"  Execution:   {exec_status}")
    stab_line = f"  PIGEAN stability: {_stability_display(stab_status)}"
    if stab_source:
        stab_line += f" ({stab_source})"
    lines.append(stab_line)
    lines.append(f"  Formal MCMC convergence: NOT FORMALLY ASSESSED")
    lines.append("")

    # Analysis settings
    lines.append("ANALYSIS SETTINGS")
    lines.append("-" * 40)
    analysis = manifest.get("analysis", {})
    lines.append(f"  Analysis type:  {analysis.get('analysis_type', 'unknown')}")
    lines.append(f"  Preset:         {analysis.get('preset', 'unknown')}")
    lines.append(f"  Genome build:   {analysis.get('genome_build', 'unknown')}")
    lines.append(f"  Gene sets:      {analysis.get('gene_sets', 'unknown')}")
    lines.append("")

    # Custom gene-set info (Phase 2)
    custom_gs = manifest.get("custom_gene_sets")
    if custom_gs:
        lines.append("CUSTOM GENE SETS")
        lines.append("-" * 40)
        lines.append(f"  Format:  {custom_gs.get('format', 'unknown')}")
        lines.append(f"  Action:  {custom_gs.get('action', 'unknown')}")
        lines.append(f"  Files:   {custom_gs.get('file_count', 0)}")
        for fname in custom_gs.get("files", {}):
            lines.append(f"    - {fname}")
        lines.append("")

    # Background info (Phase 2)
    bg = manifest.get("background")
    if bg and bg.get("provided"):
        lines.append("BACKGROUND GENES")
        lines.append("-" * 40)
        bg_qc = bg.get("qc")
        if bg_qc:
            lines.append(f"  Total genes:              {bg_qc.get('total_count', '?')}")
            lines.append(f"  Recognized in gene map:   {bg_qc.get('mapped_count', '?')}")
            lines.append(f"  Unmapped:                 {bg_qc.get('unmapped_count', '?')}")
            missing = bg_qc.get("input_genes_missing", [])
            if missing:
                lines.append(f"  Input genes missing from background: {', '.join(missing)}")
        lines.append("")

    # Evidence info (Phase 3 — non-positive-controls modes)
    evidence_info = manifest.get("evidence")
    if evidence_info:
        lines.append("EVIDENCE INPUT")
        lines.append("-" * 40)
        lines.append(f"  Evidence type:     {evidence_info.get('type', 'unknown')}")
        lines.append(f"  Original rows:     {evidence_info.get('original_rows', '?')}")
        lines.append(f"  Normalized rows:   {evidence_info.get('normalized_rows', '?')}")
        lines.append(f"  Removed rows:      {evidence_info.get('removed_rows', '?')}")
        lines.append(f"  Duplicates:        {evidence_info.get('duplicates_removed', '?')}")
        if evidence_info.get("gene_column"):
            lines.append(f"  Gene column:       {evidence_info['gene_column']}")
        if evidence_info.get("score_column"):
            lines.append(f"  Score column:      {evidence_info['score_column']}")
        transformations = evidence_info.get("transformations", [])
        if transformations:
            lines.append(f"  Transformations:")
            for t in transformations:
                lines.append(f"    - {t}")
        lines.append("")

    # Input QC summary
    analysis_type = config.get("analysis", "positive-controls")
    is_gwas = (analysis_type == "gwas")

    if is_gwas:
        lines.append("INPUT GWAS QC")
        lines.append("-" * 40)
        genome_build = config.get("genome_build", "unknown")
        lines.append(f"  Evidence type:  GWAS summary statistics (SNP-level)")
        lines.append(f"  Genome build:   {genome_build}")
        lines.append(f"  Gene QC:        Not applicable — SNP-to-gene mapping is")
        lines.append(f"                  performed by the engine using reference files.")
        lines.append("")
        if qc_result is not None and qc_result.issues:
            lines.append("  QC notes:")
            for issue in qc_result.issues:
                lines.append(f"    - {issue}")
            lines.append("")

        gwas_cols_specified = []
        for key in ["gwas_chrom_col", "gwas_pos_col", "gwas_p_col",
                     "gwas_beta_col", "gwas_se_col", "gwas_n_col",
                     "gwas_freq_col", "gwas_locus_col"]:
            val = config.get(key)
            if val:
                gwas_cols_specified.append((key, val))
        if gwas_cols_specified:
            lines.append("  Column overrides:")
            for key, val in gwas_cols_specified:
                flag = key.replace("_", "-")
                lines.append(f"    --{flag} = {val}")
            lines.append("")
        gwas_n = config.get("gwas_n")
        if gwas_n is not None:
            lines.append(f"  Global sample size (--gwas-n): {gwas_n}")
            lines.append("")

    elif hasattr(qc_result, "input_count"):
        lines.append("INPUT GENE QC")
        lines.append("-" * 40)
        genome_build = config.get("genome_build", "hg19")
        lines.append(f"  Input genes:              {qc_result.input_count}")
        lines.append(f"  Recognized in gene map:   {qc_result.recognized_count}")
        lines.append(f"  Unresolved:               {qc_result.unresolved_count}")
        lines.append(f"  With coordinates ({genome_build}):  {qc_result.with_coordinates_count}")
        lines.append(f"  In annotation gene sets:  {qc_result.with_annotations_count}")
        lines.append("")

        if qc_result.per_gene:
            lines.append("  Per-gene status:")
            for g in qc_result.per_gene:
                status = "recognized" if g["recognized"] else "UNRESOLVED"
                annot = g.get("annotation_memberships", 0)
                lines.append(f"    {g['gene']:20s}  {status:12s}  annotations={annot}")
            lines.append("")

        if qc_result.issues:
            lines.append("  QC warnings/issues:")
            for issue in qc_result.issues:
                lines.append(f"    - {issue}")
            lines.append("")

    elif qc_result is not None and hasattr(qc_result, "total_genes"):
        lines.append("INPUT EVIDENCE QC")
        lines.append("-" * 40)
        lines.append(f"  Evidence type:            {getattr(qc_result, 'evidence_type', 'unknown')}")
        lines.append(f"  Total genes:              {qc_result.total_genes}")
        lines.append(f"  Recognized in gene map:   {qc_result.recognized_genes}")
        lines.append(f"  Unresolved:               {qc_result.unresolved_genes}")
        if qc_result.duplicate_genes:
            lines.append(f"  Duplicate genes:          {qc_result.duplicate_genes}")
        if qc_result.missing_values:
            lines.append(f"  Missing values:           {qc_result.missing_values}")
        if qc_result.invalid_values:
            lines.append(f"  Invalid values:           {qc_result.invalid_values}")
        lines.append("")

        if qc_result.issues:
            lines.append("  QC warnings/issues:")
            for issue in qc_result.issues:
                lines.append(f"    - {issue}")
            lines.append("")
    else:
        lines.append("INPUT QC")
        lines.append("-" * 40)
        lines.append("  (QC not performed — input validation failed)")
        lines.append("")

    # Validation issues
    val_issues = manifest.get("validation", {}).get("issues", [])
    if val_issues:
        lines.append("VALIDATION ISSUES")
        lines.append("-" * 40)
        for issue in val_issues:
            lines.append(f"  - {issue}")
        lines.append("")

    # Execution details
    if exec_status == "COMPLETED":
        lines.append("EXECUTION")
        lines.append("-" * 40)
        lines.append(f"  Engine exit code: {manifest.get('execution', {}).get('exit_code', '?')}")
        duration = manifest.get("timing", {}).get("duration_seconds")
        if duration is not None:
            lines.append(f"  Duration: {duration:.1f} seconds")
        lines.append("")

        # Output files
        lines.append("OUTPUT FILES")
        lines.append("-" * 40)
        for fname in ["gs.out", "gss.out", "ggss.out", "p.out"]:
            fpath = os.path.join(output_dir, fname)
            if os.path.exists(fpath):
                with open(fpath) as f:
                    row_count = sum(1 for _ in f) - 1
                lines.append(f"  {fname:20s}  {row_count:>8d} data rows")
            else:
                lines.append(f"  {fname:20s}  MISSING")
        lines.append("")

        # ── Phase 5: PIGEAN Stability Assessment ──
        if convergence_result is not None:
            _add_stability_section(lines, convergence_result)

        # ── Phase 5: Interpretation sections ──
        if interpretation_result is not None:
            _add_gene_summary_section(lines, interpretation_result)
            _add_gene_set_summary_section(lines, interpretation_result)
            _add_gene_pathway_section(lines, interpretation_result)
            _add_metric_guide_section(lines, interpretation_result)
            _add_caveats_section(lines, interpretation_result)
        else:
            # Legacy interpretation notes when no parsed output
            _add_legacy_interpretation_notes(lines, config)

    elif exec_status == "NOT_RUN":
        lines.append("ENGINE NOT EXECUTED")
        lines.append("-" * 40)
        lines.append("  Validation failed — priors.py was not run.")
        lines.append("  Review validation issues above.")
        lines.append("")

    elif exec_status == "FAILED":
        lines.append("ENGINE FAILED")
        lines.append("-" * 40)
        exit_code = manifest.get("execution", {}).get("exit_code", "?")
        lines.append(f"  Engine exit code: {exit_code}")
        lines.append(f"  See pigean_run.log for details.")
        lines.append("")
        # Still show stability if available
        if convergence_result is not None:
            _add_stability_section(lines, convergence_result)

    # Full paths
    lines.append("FILE LOCATIONS")
    lines.append("-" * 40)
    qc_file = ("input_gene_qc.tsv" if analysis_type == "positive-controls"
                else "input_evidence_qc.tsv")
    for fname in ["resolved_config.json", qc_file,
                   "priors_command.txt", "run_manifest.json",
                   "pigean_run.log", "convergence.json",
                   "report.txt", "report.html"]:
        if fname in ("report.txt", "report.html"):
            exists = "✓"
        else:
            fpath = os.path.join(output_dir, fname)
            exists = "✓" if os.path.exists(fpath) else "✗"
        lines.append(f"  {exists} {fname}")
    lines.append("")

    lines.append("=" * 60)

    report_path = os.path.join(output_dir, "report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return report_path


# ── Display helpers ──────────────────────────────────────────────────


def _stability_display(status_str):
    """Human-friendly display for stability status strings."""
    mapping = {
        "PIGEAN_STABILITY_CRITERION_MET": "MET",
        "PIGEAN_STABILITY_CRITERION_UNCERTAIN": "UNCERTAIN",
        "PIGEAN_STABILITY_CRITERION_NOT_MET": "NOT MET",
        "NOT_ASSESSED": "NOT ASSESSED",
        "TRACE_NOT_AVAILABLE": "NOT ASSESSED (trace unavailable)",
        "ENGINE_FAILED": "NOT ASSESSED (engine failed)",
    }
    return mapping.get(status_str, status_str)


# ── Phase 5 section helpers ──────────────────────────────────────────


def _add_stability_section(lines, convergence_result):
    """Add PIGEAN STABILITY ASSESSMENT section to report."""
    lines.append("PIGEAN STABILITY ASSESSMENT")
    lines.append("-" * 40)

    stab_status = convergence_result["stability_status"].value
    source = convergence_result["assessment_source"]
    iter_status = convergence_result["iteration_status"].value
    details = convergence_result.get("details", {})

    lines.append(f"  PIGEAN stability criterion:  "
                 f"{_stability_display(stab_status)}")
    lines.append("")
    lines.append(f"  Method:")
    lines.append(f"    Maximum fractional SEM across cross-chain running-average")
    lines.append(f"    gene log-posterior-odds (engine built-in criterion)")
    lines.append("")

    sem_ratio = details.get("final_max_sem_ratio")
    threshold = details.get("max_frac_sem_threshold", 0.01)
    if sem_ratio is not None:
        lines.append(f"  Observed max fractional SEM:  {sem_ratio:.6g}")
    lines.append(f"  Threshold (--max-frac-sem):   {threshold}")
    lines.append("")

    lines.append(f"  Assessment source:  {source}")
    lines.append(f"  Iteration status:   "
                 f"{iter_status.replace('_', ' ')}")
    lines.append("")
    lines.append(f"  Formal MCMC convergence:  NOT FORMALLY ASSESSED")
    lines.append(f"  Reason:")
    lines.append(f"    The frozen PIGEAN engine does not retain sufficient")
    lines.append(f"    independent chain traces for standard multi-chain")
    lines.append(f"    diagnostics such as R-hat or effective sample size.")
    lines.append("")

    evidence = convergence_result.get("evidence", [])
    if evidence:
        lines.append("  Evidence:")
        for e in evidence:
            lines.append(f"    - {e}")

    warnings = convergence_result.get("warnings", [])
    if warnings:
        lines.append("  Warnings:")
        for w in warnings:
            lines.append(f"    - {w}")

    lines.append("")


def _add_gene_summary_section(lines, interpretation_result):
    """Add TOP PRIORITIZED GENES section."""
    gs = interpretation_result.get("gene_summary")
    if gs is None:
        return

    lines.append("TOP PRIORITIZED GENES (by model-derived posterior probability)")
    lines.append("-" * 40)
    lines.append("  Note: combined_D is the model-derived posterior probability of")
    lines.append("  gene-disease association: exp(combined)/(1+exp(combined)).")
    lines.append("  This is NOT a validated probability of disease causation.")
    lines.append("")

    # Input genes
    input_genes = gs.get("input_genes", [])
    if input_genes:
        lines.append("  INPUT GENES (positive controls — high scores expected by design):")
        _add_gene_table(lines, input_genes)
        lines.append("")

    # Candidate genes
    candidate_genes = gs.get("candidate_genes", [])
    if candidate_genes:
        lines.append("  CANDIDATE GENES (model-prioritized, not input):")
        _add_gene_table(lines, candidate_genes)
        lines.append("")

    if not input_genes and not candidate_genes:
        lines.append("  No genes to display.")
        lines.append("")


def _add_gene_table(lines, genes):
    """Add a formatted gene table."""
    header = (f"  {'Rank':>4s}  {'Gene':20s}  {'combined_D':>10s}  "
              f"{'prior':>8s}  {'log_bf':>8s}  {'N':>5s}  Location")
    lines.append(header)
    for g in genes:
        loc = g.get("location", "?")
        lines.append(
            f"  {g['rank']:4d}  {g['gene']:20s}  "
            f"{g['combined_D']:10.4f}  {g['prior']:8.3f}  "
            f"{g['log_bf']:8.3f}  {g['N']:5d}  {loc}"
        )


def _add_gene_set_summary_section(lines, interpretation_result):
    """Add TOP ENRICHED GENE SETS section."""
    gss = interpretation_result.get("gene_set_summary")
    if gss is None:
        return

    lines.append("TOP ENRICHED GENE SETS (by posterior mean effect)")
    lines.append("-" * 40)
    lines.append("  Note: beta is the posterior mean effect of gene set on gene priors,")
    lines.append("  corrected for LD between gene sets. avg_postp is the posterior")
    lines.append("  inclusion probability (0-1).")
    lines.append("")

    top_gs = gss.get("top_gene_sets", [])
    if top_gs:
        header = (f"  {'Rank':>4s}  {'Gene_Set':44s}  {'N':>5s}  "
                  f"{'beta':>10s}  {'avg_postp':>10s}  {'P':>10s}")
        lines.append(header)
        for g in top_gs:
            gs_name = g["gene_set"]
            if len(gs_name) > 44:
                gs_name = gs_name[:41] + "..."
            lines.append(
                f"  {g['rank']:4d}  {gs_name:44s}  {g['N']:5d}  "
                f"{g['beta']:10.5f}  {g['avg_postp']:10.5f}  "
                f"{g['P']:10.4g}"
            )
    lines.append("")


def _add_gene_pathway_section(lines, interpretation_result):
    """Add GENE-PATHWAY LINKS section."""
    links = interpretation_result.get("gene_pathway_links")
    if links is None:
        return

    by_gene = links.get("by_gene", {})
    if not by_gene:
        return

    lines.append("GENE-PATHWAY LINKS (top genes and their contributing gene sets)")
    lines.append("-" * 40)

    for gene_name, info in by_gene.items():
        status = info.get("status", "CANDIDATE")
        combined_D = info.get("combined_D", 0)
        lines.append(f"  {gene_name} (combined_D={combined_D:.4f}, {status}):")
        for gs_link in info.get("contributing_gene_sets", []):
            gs_name = gs_link["gene_set"]
            if len(gs_name) > 44:
                gs_name = gs_name[:41] + "..."
            lines.append(f"    {gs_name:44s}  beta={gs_link['beta']:.5f}")
        lines.append("")


def _add_metric_guide_section(lines, interpretation_result):
    """Add METRIC INTERPRETATION GUIDE section."""
    contracts = interpretation_result.get("interpretation_contracts", {})
    if not contracts:
        return

    lines.append("METRIC INTERPRETATION GUIDE")
    lines.append("-" * 40)

    guide = [
        ("combined_D", "Model-derived posterior probability [0,1]. Higher = stronger"
         " model evidence for gene-disease association. Computed as"
         " exp(combined)/(1+exp(combined)). NOT a validated probability."),
        ("prior", "Gene-level prior log-odds from gene-set model. Higher ="
         " more annotation evidence."),
        ("log_bf", "Gene-level log Bayes factor from input data. Higher ="
         " stronger data evidence."),
        ("combined", "prior + log_bf (log-odds scale). Combined evidence."),
        ("beta (gss)", "Posterior mean effect of gene set, LD-corrected."
         " Higher |value| = stronger."),
        ("avg_postp", "Posterior inclusion probability for gene set [0,1]."
         " Higher = more likely non-zero."),
        ("P (gss)", "Marginal regression p-value (univariate, initial filter)."
         " May disagree with Gibbs beta."),
    ]
    for name, desc in guide:
        lines.append(f"  {name:16s} {desc}")

    lines.append("")


def _add_caveats_section(lines, interpretation_result):
    """Add CAVEATS section."""
    caveats = interpretation_result.get("caveats", [])
    if not caveats:
        return

    lines.append("CAVEATS")
    lines.append("-" * 40)
    for caveat in caveats:
        lines.append(f"  - {caveat}")

    # Always add the stability/convergence distinction caveat
    lines.append("")
    lines.append("  NOTE ON STABILITY vs. CONVERGENCE:")
    lines.append("  PIGEAN's engine-specific stability criterion (max fractional SEM)")
    lines.append("  measures cross-chain agreement of gene posterior estimates. This is")
    lines.append("  the engine's built-in diagnostic from the original PIGEAN codebase.")
    lines.append("  Formal multi-chain MCMC convergence (R-hat / effective sample size)")
    lines.append("  was not assessed because the frozen engine does not retain the")
    lines.append("  independent chain traces needed for such diagnostics.")
    lines.append("")


def _add_legacy_interpretation_notes(lines, config):
    """Legacy interpretation notes for runs without Phase 5 parsing."""
    analysis_type = config.get("analysis", "positive-controls")

    lines.append("INTERPRETATION NOTES")
    lines.append("-" * 40)
    if analysis_type == "positive-controls":
        lines.append("  - Positive control genes are INPUT genes, not independently")
        lines.append("    discovered genes. Their high priors are expected by design.")
        lines.append("  - Additional genes with elevated priors are prioritized")
        lines.append("    candidates based on gene-set co-membership patterns.")
    elif analysis_type == "gene-bayes-factor":
        lines.append("  - Input evidence consists of natural-log Bayes factors.")
    elif analysis_type == "gene-z-score":
        lines.append("  - Input values are treated as log-odds by the engine.")
    elif analysis_type == "gene-percentile":
        lines.append("  - Input values are used to rank genes.")
    elif analysis_type == "exome":
        lines.append("  - Input is gene-level exome/rare-variant association statistics.")
    elif analysis_type == "gwas":
        genome_build = config.get("genome_build", "unknown")
        lines.append("  - Input is GWAS summary statistics (SNP-level evidence).")
        lines.append(f"    Genome build: {genome_build}.")
    lines.append("  - Stability has NOT been assessed for this run.")
    lines.append("  - Gene priors are stochastic; re-running will produce")
    lines.append("    numerically different but correlated results.")
    lines.append("  - See run_manifest.json for full provenance.")
    lines.append("")


# ── HTML Report ──────────────────────────────────────────────────────


def generate_report_html(config, manifest, output_dir,
                         interpretation_result=None,
                         convergence_result=None):
    """Generate a self-contained HTML report alongside the text report.

    Uses inline CSS only (no external dependencies).
    """
    analysis_type = config.get("analysis", "positive-controls")
    exec_status = manifest.get("execution", {}).get("status", "UNKNOWN")

    # Stability badge
    if convergence_result is not None:
        stab_status = convergence_result["stability_status"].value
        stab_source = convergence_result["assessment_source"]
    else:
        stab_status = "NOT_ASSESSED"
        stab_source = None

    badge_colors = {
        "PIGEAN_STABILITY_CRITERION_MET": ("#2d6a2d", "#e8f5e8"),
        "PIGEAN_STABILITY_CRITERION_UNCERTAIN": ("#8a6d00", "#fff8e0"),
        "PIGEAN_STABILITY_CRITERION_NOT_MET": ("#8b1a1a", "#fde8e8"),
        "NOT_ASSESSED": ("#666", "#f0f0f0"),
        "TRACE_NOT_AVAILABLE": ("#666", "#f0f0f0"),
        "ENGINE_FAILED": ("#8b1a1a", "#fde8e8"),
    }
    badge_fg, badge_bg = badge_colors.get(stab_status, ("#666", "#f0f0f0"))
    stab_display = _stability_display(stab_status)

    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html lang='en'><head><meta charset='UTF-8'>")
    parts.append("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    parts.append("<title>PIGEAN Run Report</title>")
    parts.append("<style>")
    parts.append(_get_report_css())
    parts.append("</style></head><body>")
    parts.append("<div class='container'>")

    # Header
    parts.append("<h1>PIGEAN Run Report</h1>")
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    parts.append(f"<p class='meta'>Generated: {ts} | "
                 f"Wrapper: {manifest.get('wrapper_version', 'unknown')}</p>")

    # Status badges
    parts.append("<div class='status-row'>")
    parts.append(f"<span class='badge' style='background:{badge_bg};color:{badge_fg}'>"
                 f"PIGEAN Stability: {stab_display}</span>")
    parts.append(f"<span class='badge'>Execution: {exec_status}</span>")
    parts.append(f"<span class='badge'>Analysis: {analysis_type}</span>")
    parts.append(f"<span class='badge' style='background:#f0f0f0;color:#666'>"
                 f"Formal MCMC Convergence: NOT FORMALLY ASSESSED</span>")
    parts.append("</div>")

    # Stability details
    if convergence_result is not None:
        parts.append("<h2>PIGEAN Stability Assessment</h2>")
        parts.append(f"<p><strong>PIGEAN stability criterion:</strong> "
                     f"{stab_display}</p>")
        parts.append(f"<p><strong>Assessment source:</strong> {stab_source}</p>")

        details = convergence_result.get("details", {})
        sem_ratio = details.get("final_max_sem_ratio")
        threshold = details.get("max_frac_sem_threshold", 0.01)

        parts.append("<table class='small'>")
        parts.append("<tr><th>Metric</th><th>Value</th></tr>")
        if sem_ratio is not None:
            parts.append(f"<tr><td>Max fractional SEM</td>"
                         f"<td>{sem_ratio:.6g}</td></tr>")
        parts.append(f"<tr><td>Threshold (--max-frac-sem)</td>"
                     f"<td>{threshold}</td></tr>")
        iter_status = convergence_result["iteration_status"].value
        parts.append(f"<tr><td>Iteration status</td>"
                     f"<td>{iter_status.replace('_', ' ')}</td></tr>")
        parts.append("</table>")

        parts.append("<p><strong>Formal MCMC convergence:</strong> "
                     "NOT FORMALLY ASSESSED</p>")
        parts.append("<p class='note'>The frozen PIGEAN engine does not retain "
                     "sufficient independent chain traces for standard "
                     "multi-chain diagnostics such as R-hat or effective "
                     "sample size. The stability criterion above is the "
                     "engine's built-in max-fractional-SEM check, not a "
                     "formal convergence diagnostic.</p>")

        evidence = convergence_result.get("evidence", [])
        if evidence:
            parts.append("<p><strong>Evidence:</strong></p><ul>")
            for e in evidence:
                parts.append(f"<li>{_html_escape(e)}</li>")
            parts.append("</ul>")
        warnings = convergence_result.get("warnings", [])
        if warnings:
            parts.append("<p><strong>Warnings:</strong></p><ul>")
            for w in warnings:
                parts.append(f"<li>{_html_escape(w)}</li>")
            parts.append("</ul>")

    # Interpretation sections
    if interpretation_result is not None and exec_status == "COMPLETED":
        gs = interpretation_result.get("gene_summary")
        if gs:
            parts.append("<h2>Top Prioritized Genes</h2>")
            parts.append("<p class='note'>combined_D is the model-derived posterior "
                         "probability of gene-disease association: "
                         "exp(combined)/(1+exp(combined)). "
                         "<strong>NOT</strong> a validated probability of disease causation.</p>")

            input_genes = gs.get("input_genes", [])
            if input_genes:
                parts.append("<h3>Input Genes (positive controls — high scores expected)</h3>")
                parts.append(_gene_table_html(input_genes))

            candidate_genes = gs.get("candidate_genes", [])
            if candidate_genes:
                parts.append("<h3>Candidate Genes (model-prioritized)</h3>")
                parts.append(_gene_table_html(candidate_genes))

        gss = interpretation_result.get("gene_set_summary")
        if gss and gss.get("top_gene_sets"):
            parts.append("<h2>Top Enriched Gene Sets</h2>")
            parts.append("<p class='note'>beta: posterior mean effect (LD-corrected). "
                         "avg_postp: posterior inclusion probability.</p>")
            parts.append(_gene_set_table_html(gss["top_gene_sets"]))

        # Gene-pathway links
        links = interpretation_result.get("gene_pathway_links", {})
        by_gene = links.get("by_gene", {})
        if by_gene:
            parts.append("<h2>Gene-Pathway Links</h2>")
            for gene_name, info in by_gene.items():
                status = info.get("status", "CANDIDATE")
                combined_D = info.get("combined_D", 0)
                parts.append(f"<details><summary><strong>{_html_escape(gene_name)}</strong> "
                             f"(combined_D={combined_D:.4f}, {status})</summary>")
                gs_links = info.get("contributing_gene_sets", [])
                if gs_links:
                    parts.append("<table class='small'><tr><th>Gene Set</th>"
                                 "<th>beta</th></tr>")
                    for gl in gs_links:
                        parts.append(f"<tr><td>{_html_escape(gl['gene_set'])}</td>"
                                     f"<td>{gl['beta']:.5f}</td></tr>")
                    parts.append("</table>")
                parts.append("</details>")

        # Metric guide
        contracts = interpretation_result.get("interpretation_contracts", {})
        if contracts:
            parts.append("<h2>Metric Interpretation Guide</h2>")
            parts.append("<table><tr><th>Metric</th><th>Definition</th>"
                         "<th>Scale</th><th>Stochastic</th></tr>")
            for name, c in contracts.items():
                stoch = "Yes" if c.get("stochastic") else "No"
                parts.append(f"<tr><td><strong>{_html_escape(name)}</strong></td>"
                             f"<td>{_html_escape(c.get('definition', ''))}"
                             f"<br><em>{_html_escape(c.get('anti_interpretation', ''))}</em></td>"
                             f"<td>{_html_escape(c.get('scale', ''))}</td>"
                             f"<td>{stoch}</td></tr>")
            parts.append("</table>")

        # Caveats
        caveats = interpretation_result.get("caveats", [])
        if caveats:
            parts.append("<h2>Caveats</h2><ul>")
            for c in caveats:
                parts.append(f"<li>{_html_escape(c)}</li>")
            parts.append("</ul>")

    parts.append("</div></body></html>")

    html_path = os.path.join(output_dir, "report.html")
    with open(html_path, "w") as f:
        f.write("\n".join(parts))

    return html_path


def _gene_table_html(genes):
    """Generate an HTML table for a list of gene entries."""
    rows = ["<table>",
            "<tr><th>Rank</th><th>Gene</th><th>combined_D</th>"
            "<th>prior</th><th>log_bf</th><th>N</th><th>Location</th></tr>"]
    for g in genes:
        rows.append(
            f"<tr><td>{g['rank']}</td><td>{_html_escape(g['gene'])}</td>"
            f"<td>{g['combined_D']:.4f}</td><td>{g['prior']:.3f}</td>"
            f"<td>{g['log_bf']:.3f}</td><td>{g['N']}</td>"
            f"<td>{_html_escape(g.get('location', '?'))}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _gene_set_table_html(gene_sets):
    """Generate an HTML table for gene set entries."""
    rows = ["<table>",
            "<tr><th>Rank</th><th>Gene Set</th><th>N</th>"
            "<th>beta</th><th>avg_postp</th><th>P</th></tr>"]
    for g in gene_sets:
        rows.append(
            f"<tr><td>{g['rank']}</td>"
            f"<td>{_html_escape(g['gene_set'])}</td>"
            f"<td>{g['N']}</td><td>{g['beta']:.5f}</td>"
            f"<td>{g['avg_postp']:.5f}</td><td>{g['P']:.4g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _html_escape(s):
    """Simple HTML escaping."""
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _get_report_css():
    """Inline CSS for the HTML report."""
    return """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       line-height: 1.6; color: #333; margin: 0; padding: 20px; background: #f8f9fa; }
.container { max-width: 1100px; margin: 0 auto; background: #fff;
             padding: 30px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
h1 { border-bottom: 2px solid #2c5282; padding-bottom: 10px; color: #1a365d; }
h2 { color: #2c5282; margin-top: 30px; border-bottom: 1px solid #e2e8f0; padding-bottom: 5px; }
h3 { color: #4a5568; }
.meta { color: #718096; font-size: 0.9em; }
.status-row { display: flex; gap: 10px; margin: 15px 0; flex-wrap: wrap; }
.badge { padding: 4px 12px; border-radius: 4px; font-size: 0.85em;
         background: #edf2f7; color: #4a5568; font-weight: 500; }
.note { background: #fffbeb; border-left: 3px solid #d69e2e; padding: 8px 12px;
        font-size: 0.9em; margin: 10px 0; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 0.85em; }
th, td { padding: 6px 10px; text-align: left; border-bottom: 1px solid #e2e8f0; }
th { background: #f7fafc; font-weight: 600; color: #2d3748; }
tr:hover { background: #f7fafc; }
table.small { width: auto; }
details { margin: 5px 0; }
summary { cursor: pointer; padding: 4px; }
summary:hover { background: #f7fafc; }
ul { padding-left: 20px; }
em { color: #718096; font-size: 0.9em; }
"""
