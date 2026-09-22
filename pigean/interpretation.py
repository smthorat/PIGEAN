"""
pigean/interpretation.py — Scientific interpretation of PIGEAN outputs.

Phase 5: Generates structured interpretation with strict language constraints.
Only VERIFIED metrics are interpreted. Unverified metrics are not displayed
as interpreted results.

VERIFIED metrics: combined_D, prior, log_bf, combined, N, Gene, Chrom/Start/End,
                  positive_control, beta (gss), avg_postp (gss), P (gss),
                  Gene_Set, label, N (gss), beta/weight (ggss)

Language constraints enforced:
  NEVER: "causal", "disease-causing", "validated", "confirmed",
         "true probability", "newly discovered", "independently identified"
  ALWAYS: "prioritized", "model-derived", "candidate", "associated", "enriched"
"""


# ─── Forbidden language terms ──────────────────────────────────────────

FORBIDDEN_TERMS = [
    "causal",
    "disease-causing",
    "validated",
    "confirmed",
    "true probability",
    "newly discovered",
    "independently identified",
]


# ─── Interpretation contracts ──────────────────────────────────────────

INTERPRETATION_CONTRACTS = {
    "combined_D": {
        "definition": (
            "Model-derived posterior probability of gene-disease association"
        ),
        "formula": (
            "exp(combined) / (1 + exp(combined)) "
            "where combined = prior + log_bf"
        ),
        "scale": "[0, 1]",
        "direction": (
            "Higher values indicate stronger model evidence for association"
        ),
        "anti_interpretation": (
            "This is NOT a measure of whether a gene truly contributes to "
            "disease. It is the model's assessment given gene-set annotations "
            "and input data."
        ),
        "stochastic": True,
    },
    "prior": {
        "definition": (
            "Gene-level prior log-odds from gene-set annotation model"
        ),
        "formula": (
            "X * beta / scale_factors, averaged across Gibbs chains "
            "and iterations"
        ),
        "scale": "Log-odds (unbounded)",
        "direction": (
            "Higher values indicate more annotation evidence from gene sets"
        ),
        "anti_interpretation": (
            "Does not incorporate observed data (that is log_bf). Reflects "
            "gene-set membership patterns only."
        ),
        "stochastic": True,
    },
    "log_bf": {
        "definition": "Gene-level log Bayes factor from observed data",
        "formula": (
            "Computed from input evidence, method depends on analysis type"
        ),
        "scale": "Log-odds (unbounded)",
        "direction": "Higher values indicate stronger data evidence",
        "anti_interpretation": (
            "In positive-controls mode, this IS the input signal itself, "
            "not independent evidence."
        ),
        "stochastic": False,
    },
    "combined": {
        "definition": "Combined evidence score: prior + log_bf",
        "formula": "prior + log_bf",
        "scale": "Log-odds (unbounded)",
        "direction": "Higher values indicate stronger combined evidence",
        "anti_interpretation": (
            "Sum of annotation-based prior and data-based log Bayes factor."
        ),
        "stochastic": True,
    },
    "beta_gss": {
        "definition": (
            "Posterior mean effect of gene set on gene priors, corrected "
            "for LD between gene sets"
        ),
        "formula": (
            "Average across Gibbs chains of posterior mean beta "
            "(spike-and-slab), divided by scale factor"
        ),
        "scale": "Per-gene log-odds contribution (unbounded)",
        "direction": (
            "Higher absolute value indicates stronger gene-set effect"
        ),
        "anti_interpretation": (
            "Not a univariate measure. Conditioned on all other gene sets "
            "via V matrix correction."
        ),
        "stochastic": True,
    },
    "avg_postp": {
        "definition": "Posterior inclusion probability for gene set",
        "formula": (
            "Average posterior probability of non-zero effect across "
            "Gibbs iterations"
        ),
        "scale": "[0, 1]",
        "direction": (
            "Higher values indicate the gene set more likely has a "
            "non-zero effect"
        ),
        "anti_interpretation": (
            "Not a p-value. Not directly comparable to the P column."
        ),
        "stochastic": True,
    },
    "P_gss": {
        "definition": (
            "Marginal regression p-value from initial logistic filter"
        ),
        "formula": "2 * norm.cdf(-|Z|) where Z = beta_tilde / SE",
        "scale": "[0, 1]",
        "direction": (
            "Lower values indicate more statistical significance in "
            "univariate analysis"
        ),
        "anti_interpretation": (
            "From the INITIAL filter, not the Gibbs sampler. May disagree "
            "with multivariate beta."
        ),
        "stochastic": True,
    },
}


# ─── Public API ────────────────────────────────────────────────────────

def interpret_results(parsed_gs, parsed_gss, parsed_ggss, parsed_params,
                      convergence_result, analysis_type,
                      top_n_genes=20, top_n_gene_sets=20):
    """Generate structured interpretation from parsed outputs.

    Only VERIFIED metrics are interpreted.  Unverified metrics are omitted.

    Args:
        parsed_gs: Output of parse_gene_stats(), or None.
        parsed_gss: Output of parse_gene_set_stats(), or None.
        parsed_ggss: Output of parse_gene_gene_set_stats(), or None.
        parsed_params: Output of parse_params(), or None.
        convergence_result: Dict from assess_convergence(), or None.
        analysis_type: Analysis mode string (e.g. "positive-controls").
        top_n_genes: Max genes per group (INPUT / CANDIDATE).
        top_n_gene_sets: Max gene sets to include.

    Returns:
        dict  –  structured interpretation (see module docstring).
    """
    result = {
        "gene_summary": None,
        "gene_set_summary": None,
        "gene_pathway_links": None,
        "convergence_note": _convergence_note(convergence_result),
        "interpretation_contracts": INTERPRETATION_CONTRACTS,
        "caveats": [],
        "analysis_type": analysis_type,
    }

    has_positive_control = False

    # ── Gene summary ───────────────────────────────────────────────
    if parsed_gs is not None:
        has_positive_control = parsed_gs.get("has_positive_control", False)
        genes = parsed_gs.get("genes", [])

        input_genes, candidate_genes = _separate_input_and_candidate_genes(
            genes, has_positive_control,
        )

        result["gene_summary"] = {
            "input_genes": [
                _format_gene_entry(g, rank=i + 1, status="INPUT")
                for i, g in enumerate(input_genes[:top_n_genes])
            ],
            "candidate_genes": [
                _format_gene_entry(g, rank=i + 1, status="CANDIDATE")
                for i, g in enumerate(candidate_genes[:top_n_genes])
            ],
            "total_genes": parsed_gs.get("row_count", len(genes)),
            "note": _gene_summary_note(analysis_type, has_positive_control),
        }

    # ── Gene-set summary ───────────────────────────────────────────
    if parsed_gss is not None:
        gene_sets = parsed_gss.get("gene_sets", [])
        result["gene_set_summary"] = {
            "top_gene_sets": [
                _format_gene_set_entry(gs, rank=i + 1)
                for i, gs in enumerate(gene_sets[:top_n_gene_sets])
            ],
            "total_gene_sets": parsed_gss.get("row_count", len(gene_sets)),
            "note": _gene_set_summary_note(),
        }

    # ── Gene ↔ pathway links ──────────────────────────────────────
    if parsed_ggss is not None and parsed_gs is not None and parsed_gss is not None:
        # Collect names of top genes (both groups) and top gene sets
        top_gene_names = set()
        if result["gene_summary"] is not None:
            for g in result["gene_summary"]["input_genes"]:
                top_gene_names.add(g["gene"])
            for g in result["gene_summary"]["candidate_genes"]:
                top_gene_names.add(g["gene"])

        top_gs_names = set()
        if result["gene_set_summary"] is not None:
            for gs in result["gene_set_summary"]["top_gene_sets"]:
                top_gs_names.add(gs["gene_set"])

        result["gene_pathway_links"] = _build_gene_pathway_links(
            parsed_ggss, top_gene_names, top_gs_names,
            parsed_gs, parsed_gss,
        )

    # ── Caveats ────────────────────────────────────────────────────
    result["caveats"] = _generate_caveats(
        analysis_type, convergence_result, has_positive_control,
    )

    return result


# ─── Gene helpers ──────────────────────────────────────────────────────

def _separate_input_and_candidate_genes(genes, has_positive_control):
    """Split genes into input (positive_control > 0) and candidate groups.

    Both lists are pre-sorted by combined_D descending (parser does this).
    If there is no positive_control column, all genes are candidates.

    Returns:
        (input_genes, candidate_genes)  — two lists of gene dicts.
    """
    if not has_positive_control:
        return [], genes

    input_genes = []
    candidate_genes = []
    for g in genes:
        pc = g.get("positive_control")
        if pc is not None and _is_positive(pc):
            input_genes.append(g)
        else:
            candidate_genes.append(g)
    return input_genes, candidate_genes


def _is_positive(value):
    """Return True when the positive_control field signals an input gene."""
    try:
        return float(value) > 0
    except (TypeError, ValueError):
        return False


def _format_gene_entry(gene_dict, rank, status):
    """Extract only VERIFIED fields from a gene dict."""
    chrom = gene_dict.get("Chrom", "")
    start = gene_dict.get("Start", "")
    end = gene_dict.get("End", "")
    location = f"{chrom}:{start}-{end}" if chrom else ""

    return {
        "rank": rank,
        "gene": gene_dict.get("Gene", ""),
        "combined_D": _safe_float(gene_dict.get("combined_D")),
        "prior": _safe_float(gene_dict.get("prior")),
        "log_bf": _safe_float(gene_dict.get("log_bf")),
        "combined": _safe_float(gene_dict.get("combined")),
        "N": _safe_int(gene_dict.get("N")),
        "location": location,
        "status": status,
    }


def _gene_summary_note(analysis_type, has_positive_control):
    """Interpretation note for the gene summary section."""
    note = (
        "combined_D represents the model-derived posterior probability of "
        "gene-disease association, integrating gene-set annotations (prior) "
        "and observed data (log_bf). Higher values indicate stronger model "
        "evidence. This is NOT an externally verified probability."
    )
    if has_positive_control:
        note += (
            " INPUT genes were provided as positive controls; their elevated "
            "scores are expected by design and should not be interpreted as "
            "independent prioritization."
        )
    if analysis_type == "gwas":
        note += (
            " For GWAS input, log_bf is derived from SNP-to-gene mapping "
            "and reflects the strength of local association signal."
        )
    return note


# ─── Gene-set helpers ──────────────────────────────────────────────────

def _format_gene_set_entry(gs_dict, rank):
    """Extract only VERIFIED fields from a gene-set dict."""
    return {
        "rank": rank,
        "gene_set": gs_dict.get("Gene_Set", ""),
        "label": gs_dict.get("label", ""),
        "N": _safe_int(gs_dict.get("N")),
        "beta": _safe_float(gs_dict.get("beta")),
        "avg_postp": _safe_float(gs_dict.get("avg_postp")),
        "P": _safe_float(gs_dict.get("P")),
    }


def _gene_set_summary_note():
    """Interpretation note for the gene-set summary section."""
    return (
        "Gene sets are ranked by the absolute value of beta, the posterior "
        "mean effect on gene priors corrected for correlation between gene "
        "sets. avg_postp is the posterior inclusion probability (0-1). P is "
        "the marginal regression p-value from the initial filter and may "
        "disagree with the multivariate (Gibbs) results."
    )


# ─── Gene ↔ pathway links ─────────────────────────────────────────────

def _build_gene_pathway_links(parsed_ggss, top_gene_names, top_gs_names,
                               parsed_gs, parsed_gss, max_links=10):
    """Build bidirectional gene ↔ gene-set links from ggss.out.

    Only builds links for genes/gene-sets already in the top-N lists.
    """
    by_gene_index = parsed_ggss.get("by_gene", {})
    by_gs_index = parsed_ggss.get("by_gene_set", {})

    # Quick lookup for combined_D and status by gene name
    gene_info = {}
    if parsed_gs is not None:
        has_pc = parsed_gs.get("has_positive_control", False)
        for g in parsed_gs.get("genes", []):
            name = g.get("Gene", "")
            if name in top_gene_names:
                pc = g.get("positive_control")
                is_input = has_pc and pc is not None and _is_positive(pc)
                gene_info[name] = {
                    "combined_D": _safe_float(g.get("combined_D")),
                    "log_bf": _safe_float(g.get("log_bf")),
                    "status": "INPUT" if is_input else "CANDIDATE",
                }

    # Quick lookup for beta by gene-set name
    gs_info = {}
    if parsed_gss is not None:
        for gs in parsed_gss.get("gene_sets", []):
            name = gs.get("Gene_Set", "")
            if name in top_gs_names:
                gs_info[name] = {
                    "beta": _safe_float(gs.get("beta")),
                }

    # Build by_gene links
    links_by_gene = {}
    for gene_name in top_gene_names:
        entries = by_gene_index.get(gene_name, [])
        # Sort by |beta| descending
        sorted_entries = sorted(
            entries, key=lambda e: abs(_safe_float(e.get("beta"))),
            reverse=True,
        )
        info = gene_info.get(gene_name, {})
        links_by_gene[gene_name] = {
            "combined_D": info.get("combined_D"),
            "status": info.get("status", "CANDIDATE"),
            "contributing_gene_sets": [
                {
                    "gene_set": e.get("gene_set", ""),
                    "beta": _safe_float(e.get("beta")),
                    "weight": _safe_float(e.get("weight")),
                }
                for e in sorted_entries[:max_links]
            ],
        }

    # Build by_gene_set links
    links_by_gs = {}
    for gs_name in top_gs_names:
        entries = by_gs_index.get(gs_name, [])
        # Sort by combined_D descending (look up from gene_info)
        def _gene_combined_d(entry):
            g = entry.get("Gene", "")
            return gene_info.get(g, {}).get("combined_D") or 0.0
        sorted_entries = sorted(entries, key=_gene_combined_d, reverse=True)
        info = gs_info.get(gs_name, {})
        links_by_gs[gs_name] = {
            "beta": info.get("beta"),
            "top_genes": [
                {
                    "gene": e.get("Gene", ""),
                    "combined_D": _safe_float(e.get("combined")),
                    "log_bf": _safe_float(e.get("log_bf")),
                }
                for e in sorted_entries[:max_links]
            ],
        }

    return {
        "by_gene": links_by_gene,
        "by_gene_set": links_by_gs,
    }


# ─── Convergence note ─────────────────────────────────────────────────

def _convergence_note(convergence_result):
    """One-line summary of stability/convergence status for interpretation."""
    if convergence_result is None:
        return "Stability was not assessed."

    status = convergence_result.get("stability_status",
                                    convergence_result.get("status",
                                                           "NOT_ASSESSED"))
    status_str = status.value if hasattr(status, "value") else str(status)
    source = convergence_result.get("assessment_source",
                                    convergence_result.get("confidence",
                                                           "INSUFFICIENT"))

    if status_str == "PIGEAN_STABILITY_CRITERION_MET":
        return (
            f"PIGEAN stability criterion: MET ({source}). "
            "Results are suitable for interpretation. "
            "Formal MCMC convergence was not assessed."
        )
    elif status_str == "PIGEAN_STABILITY_CRITERION_UNCERTAIN":
        return (
            f"PIGEAN stability criterion: UNCERTAIN ({source}). "
            "Results should be interpreted with additional care. "
            "Formal MCMC convergence was not assessed."
        )
    elif status_str == "PIGEAN_STABILITY_CRITERION_NOT_MET":
        return (
            f"PIGEAN stability criterion: NOT MET ({source}). "
            "Results may be unreliable. Consider re-running with more "
            "iterations. Formal MCMC convergence was not assessed."
        )
    elif status_str == "ENGINE_FAILED":
        return "Engine did not complete. No results to interpret."
    else:
        return f"Stability status: {status_str} ({source})."


# ─── Caveats ───────────────────────────────────────────────────────────

def _generate_caveats(analysis_type, convergence_result, has_positive_control):
    """Generate analysis-appropriate caveats."""
    caveats = [
        (
            "All stochastic metrics (combined_D, prior, beta, avg_postp) "
            "will differ across runs due to MCMC sampling. The engine does "
            "not set a fixed random seed."
        ),
        (
            "Gene priors reflect gene-set annotation patterns, not direct "
            "evidence of biological mechanism. A high combined_D means the "
            "gene is prioritized by the MODEL, not that it is established "
            "as disease-relevant."
        ),
    ]

    if has_positive_control:
        caveats.append(
            "Input genes (positive controls) are expected to have high "
            "scores because they contribute to the training signal. They "
            "should NOT be interpreted as independently prioritized."
        )

    if analysis_type == "gwas":
        caveats.append(
            "GWAS-based gene scores depend on the SNP-to-gene mapping "
            "performed by the engine. Genes near strong association signals "
            "will have higher log_bf values. Results are sensitive to the "
            "genome build and reference gene-location files used."
        )

    if convergence_result is not None:
        status = convergence_result.get("stability_status",
                                        convergence_result.get("status",
                                                               "NOT_ASSESSED"))
        status_str = status.value if hasattr(status, "value") else str(status)
        if status_str in ("PIGEAN_STABILITY_CRITERION_UNCERTAIN",
                          "PIGEAN_STABILITY_CRITERION_NOT_MET"):
            caveats.append(
                "PIGEAN's engine-specific stability criterion was not "
                "fully met. Consider re-running with additional iterations "
                "(--max-num-iter) or enabling trace output "
                "(--enable-convergence-trace) for detailed diagnostics."
            )
        elif status_str in ("NOT_ASSESSED", "TRACE_NOT_AVAILABLE"):
            caveats.append(
                "Stability was not fully assessed. For detailed "
                "diagnostics, enable trace output with "
                "--enable-convergence-trace."
            )

    return caveats


# ─── Numeric helpers ───────────────────────────────────────────────────

def _safe_float(value):
    """Convert value to float, returning None on failure."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value):
    """Convert value to int, returning None on failure."""
    if value is None:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None
