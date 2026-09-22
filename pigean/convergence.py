"""
pigean/convergence.py — Stability assessment for PIGEAN Gibbs sampler.

Phase 5: Two-tier stability diagnostics based on PIGEAN's engine-specific
max-fractional-SEM criterion.

Tier 1: From p.out parameters only (SUMMARY_BASED).
Tier 2: From p.out + log file parsing (TLOG_BASED).

IMPORTANT: COMPLETED != stability criterion met.
Execution success does NOT imply stability or convergence.

The PIGEAN engine's built-in stability criterion (--max-frac-sem) measures
the maximum normalized standard error of cross-chain running-average gene
log-posterior-odds.  This is an engine-specific stability measure, NOT
equivalent to formal multi-chain MCMC convergence diagnostics (R-hat / ESS).
"""

import json
import os
import re
from enum import Enum


class StabilityStatus(Enum):
    """PIGEAN engine-specific stability status.

    These statuses reflect whether the engine's built-in max-fractional-SEM
    criterion was satisfied.  They do NOT imply formal MCMC convergence.
    """
    PIGEAN_STABILITY_CRITERION_MET = "PIGEAN_STABILITY_CRITERION_MET"
    PIGEAN_STABILITY_CRITERION_UNCERTAIN = "PIGEAN_STABILITY_CRITERION_UNCERTAIN"
    PIGEAN_STABILITY_CRITERION_NOT_MET = "PIGEAN_STABILITY_CRITERION_NOT_MET"
    NOT_ASSESSED = "NOT_ASSESSED"
    TRACE_NOT_AVAILABLE = "TRACE_NOT_AVAILABLE"
    ENGINE_FAILED = "ENGINE_FAILED"


class FormalConvergenceStatus(Enum):
    """Formal MCMC convergence status (R-hat / ESS).

    The frozen PIGEAN engine does not retain sufficient independent chain
    traces for standard multi-chain convergence diagnostics.
    """
    NOT_FORMALLY_ASSESSED = "NOT_FORMALLY_ASSESSED"


class IterationStatus(Enum):
    """Whether the engine hit its iteration cap."""
    MAX_ITERATION_CAP_REACHED = "MAX_ITERATION_CAP_REACHED"
    STOPPED_BEFORE_CAP = "STOPPED_BEFORE_CAP"
    UNKNOWN = "UNKNOWN"


# Backward-compatible alias so existing imports do not break
ConvergenceStatus = StabilityStatus


# ── Log-line regexes ──────────────────────────────────────────────────────

# "Gibbs iteration 500: ref_val=8.1; max_sem=0.0491; max_ratio=0.00606"
_RE_GIBBS_SEM = re.compile(
    r"Gibbs iteration (\d+):\s*ref_val=([\d.eE+-]+);\s*"
    r"max_sem=([\d.eE+-]+);\s*max_ratio=([\d.eE+-]+)"
)

# "Iteration 11: max ind=...; max B=0.00218; max W=0.000228;
#  max R=1.334; avg R=1.002; num above=148;"
_RE_RHAT = re.compile(
    r"Iteration (\d+):.*?max R=([\d.eE+-]+);\s*avg R=([\d.eE+-]+);\s*"
    r"num above=(\d+);"
)

_PRECISION_MSG = "Desired Gibbs precision achieved; stopping sampling"

_FORMAL_CONVERGENCE_REASON = (
    "The frozen PIGEAN engine does not retain sufficient independent "
    "chain traces for formal R-hat/ESS-style convergence diagnostics."
)


# ── Helpers ───────────────────────────────────────────────────────────────

def _get_param(params, key, as_type=None):
    """Retrieve a parameter value from parsed params dict.

    Multi-version params are stored as lists; return the last (most recent)
    value.  Single-version params are stored as scalars.
    """
    val = params.get(key)
    if val is None:
        return None
    if isinstance(val, list):
        val = val[-1]  # last version
    if as_type is not None:
        try:
            val = as_type(val)
        except (ValueError, TypeError):
            return None
    return val


# ── Log parser ────────────────────────────────────────────────────────────

def parse_convergence_from_log(log_path):
    """Extract stability diagnostics from pigean_run.log.

    Returns a dict with parsed fields, or *None* if the log file is missing
    or unreadable.
    """
    if log_path is None or not os.path.isfile(log_path):
        return None

    precision_achieved = False
    last_gibbs_iter = None
    last_ref_val = None
    last_max_sem = None
    last_max_ratio = None
    r_hat_entries = []

    try:
        with open(log_path, "r") as fh:
            for line in fh:
                # SEM stability line
                m = _RE_GIBBS_SEM.search(line)
                if m:
                    last_gibbs_iter = int(m.group(1))
                    last_ref_val = float(m.group(2))
                    last_max_sem = float(m.group(3))
                    last_max_ratio = float(m.group(4))

                # R-hat line (informational)
                m = _RE_RHAT.search(line)
                if m:
                    r_hat_entries.append((
                        int(m.group(1)),
                        float(m.group(2)),
                        float(m.group(3)),
                        int(m.group(4)),
                    ))

                # Explicit precision message (rare — requires non-default option)
                if _PRECISION_MSG in line:
                    precision_achieved = True
    except OSError:
        return None

    final_max_r_hat = r_hat_entries[-1][1] if r_hat_entries else None
    final_avg_r_hat = r_hat_entries[-1][2] if r_hat_entries else None

    return {
        "precision_achieved": precision_achieved,
        "final_gibbs_iteration": last_gibbs_iter,
        "final_max_sem": last_max_sem,
        "final_max_ratio": last_max_ratio,
        "final_ref_val": last_ref_val,
        "r_hat_entries": r_hat_entries,
        "final_max_r_hat": final_max_r_hat,
        "final_avg_r_hat": final_avg_r_hat,
    }


# ── Main assessment ───────────────────────────────────────────────────────

def assess_convergence(parsed_params, log_path=None,
                       max_num_iter=500, max_frac_sem=0.01):
    """Assess PIGEAN engine stability using the max-fractional-SEM criterion.

    This assesses the engine's built-in stability criterion, NOT formal MCMC
    convergence.  Formal convergence (R-hat / ESS) is always reported as
    NOT_FORMALLY_ASSESSED because the frozen engine does not retain sufficient
    independent chain traces.

    Parameters
    ----------
    parsed_params : dict or None
        Output of ``parse_params()`` (has a ``"params"`` key), or *None*
        if p.out was missing / unparseable.
    log_path : str or None
        Path to ``pigean_run.log``.  When provided the log is parsed for
        SEM-ratio and R-hat diagnostics (Tier 2, TLOG_BASED).
    max_num_iter : int
        Engine's ``--max-num-iter`` default.  500 for the current engine.
    max_frac_sem : float
        Engine's ``--max-frac-sem`` threshold.  0.01 for the current engine.

    Returns
    -------
    dict with keys: stability_status, assessment_source, iteration_status,
                    formal_convergence_status, formal_convergence_reason,
                    evidence, warnings, details.

    Legacy keys ``status`` and ``confidence`` are also provided for
    backward compatibility with callers that reference them.
    """
    evidence = []
    warnings = []
    details = {
        "num_gibbs_iter": None,
        "num_gibbs_restarts": None,
        "hit_iteration_cap": False,
        "max_num_iter": max_num_iter,
        "num_chains": None,
        "log_parsed": False,
        "final_max_r_hat": None,
        "final_avg_r_hat": None,
        "final_max_sem_ratio": None,
        "precision_achieved_message": False,
        "max_frac_sem_threshold": max_frac_sem,
    }

    # ── Guard: no params at all ──────────────────────────────────────
    if parsed_params is None:
        return _build_result(
            stability_status=StabilityStatus.NOT_ASSESSED,
            assessment_source="INSUFFICIENT",
            iteration_status=IterationStatus.UNKNOWN,
            evidence=["p.out not available — cannot assess stability"],
            warnings=[],
            details=details,
        )

    params = parsed_params.get("params", {})
    if not params:
        return _build_result(
            stability_status=StabilityStatus.NOT_ASSESSED,
            assessment_source="INSUFFICIENT",
            iteration_status=IterationStatus.UNKNOWN,
            evidence=["p.out contained no parameters"],
            warnings=[],
            details=details,
        )

    # ── Tier 1: p.out parameters ─────────────────────────────────────
    num_iter = _get_param(params, "num_gibbs_iter", int)
    num_restarts = _get_param(params, "num_gibbs_restarts", int)
    num_chains = _get_param(params, "num_chains", int)

    details["num_gibbs_iter"] = num_iter
    details["num_gibbs_restarts"] = num_restarts
    details["num_chains"] = num_chains

    if num_iter is None:
        return _build_result(
            stability_status=StabilityStatus.NOT_ASSESSED,
            assessment_source="INSUFFICIENT",
            iteration_status=IterationStatus.UNKNOWN,
            evidence=["p.out missing num_gibbs_iter — cannot assess stability"],
            warnings=[],
            details=details,
        )

    # num_gibbs_iter is 0-indexed.  Value 499 with max_num_iter=500 means
    # the loop ran all 500 iterations (range(500) → 0..499).
    hit_cap = num_iter >= max_num_iter - 1
    had_restarts = (num_restarts is not None and num_restarts > 0)
    details["hit_iteration_cap"] = hit_cap

    iter_status = (IterationStatus.MAX_ITERATION_CAP_REACHED if hit_cap
                   else IterationStatus.STOPPED_BEFORE_CAP)

    # Human-readable iteration count (1-indexed for the report)
    iter_display = num_iter + 1

    if num_chains is not None:
        evidence.append(f"{num_chains} chains were used")

    if not hit_cap and not had_restarts:
        status = StabilityStatus.PIGEAN_STABILITY_CRITERION_MET
        assessment_source = "SUMMARY_BASED"
        evidence.append(
            f"Sampler stopped at iteration {iter_display} of "
            f"{max_num_iter} (did not hit cap)"
        )
        evidence.append(
            f"No restarts needed (num_gibbs_restarts="
            f"{num_restarts if num_restarts is not None else 0})"
        )
    elif hit_cap and not had_restarts:
        status = StabilityStatus.PIGEAN_STABILITY_CRITERION_UNCERTAIN
        assessment_source = "SUMMARY_BASED"
        evidence.append(
            f"Sampler ran all {max_num_iter} iterations (hit cap)"
        )
        evidence.append("No restarts needed")
        warnings.append(
            "Sampler hit the maximum iteration limit. The engine's "
            "stability criterion could not be verified from p.out alone. "
            "Consider enabling trace output for detailed diagnostics."
        )
    elif not hit_cap and had_restarts:
        status = StabilityStatus.PIGEAN_STABILITY_CRITERION_UNCERTAIN
        assessment_source = "SUMMARY_BASED"
        evidence.append(
            f"Sampler stopped at iteration {iter_display} of {max_num_iter}"
        )
        evidence.append(
            f"Required {num_restarts} restart(s) — indicates initial "
            f"stability difficulty"
        )
        warnings.append(
            f"The sampler required {num_restarts} restart(s) before "
            f"stabilizing. Hyperparameters may have been adjusted."
        )
    else:  # hit_cap and had_restarts
        status = StabilityStatus.PIGEAN_STABILITY_CRITERION_NOT_MET
        assessment_source = "SUMMARY_BASED"
        evidence.append(
            f"Sampler ran all {max_num_iter} iterations AND required "
            f"{num_restarts} restart(s)"
        )
        warnings.append(
            "Strong evidence of stability difficulty. Both the "
            "iteration cap was hit and restarts were needed. Results "
            "should be interpreted with caution."
        )

    # ── Tier 2: log parsing (upgrades assessment source) ────────────────
    log_info = parse_convergence_from_log(log_path)

    if log_info is not None:
        details["log_parsed"] = True
        details["precision_achieved_message"] = log_info["precision_achieved"]
        details["final_max_r_hat"] = log_info["final_max_r_hat"]
        details["final_avg_r_hat"] = log_info["final_avg_r_hat"]
        details["final_max_sem_ratio"] = log_info["final_max_ratio"]

        # Explicit precision message (rare but definitive)
        if log_info["precision_achieved"]:
            if status != StabilityStatus.PIGEAN_STABILITY_CRITERION_MET:
                status = StabilityStatus.PIGEAN_STABILITY_CRITERION_MET
            assessment_source = "TLOG_BASED"
            evidence.append(
                "Log file: engine reported 'Desired Gibbs precision "
                "achieved; stopping sampling'"
            )

        # SEM ratio check — the engine's built-in stability criterion
        if (log_info["final_max_ratio"] is not None
                and log_info["final_max_ratio"] < max_frac_sem):
            sem_msg = (
                f"Log file: final max fractional SEM "
                f"{log_info['final_max_ratio']:.6g} < "
                f"{max_frac_sem} threshold — PIGEAN stability "
                f"criterion met"
            )
            evidence.append(sem_msg)
            # Upgrade if we were at UNCERTAIN
            if status == StabilityStatus.PIGEAN_STABILITY_CRITERION_UNCERTAIN:
                status = StabilityStatus.PIGEAN_STABILITY_CRITERION_MET
                assessment_source = "TLOG_BASED"
                if hit_cap:
                    warnings.append(
                        "Sampler hit iteration cap but PIGEAN's "
                        "max-fractional-SEM stability criterion was "
                        "met at completion"
                    )
        elif log_info["final_max_ratio"] is not None:
            evidence.append(
                f"Log file: final max fractional SEM "
                f"{log_info['final_max_ratio']:.6g} >= "
                f"{max_frac_sem} threshold — PIGEAN stability "
                f"criterion NOT met"
            )

        # R-hat diagnostics (informational only — burn-in R-hat, not formal)
        if log_info["final_avg_r_hat"] is not None:
            avg_r = log_info["final_avg_r_hat"]
            if avg_r < 1.05:
                evidence.append(
                    f"Log file: final avg R-hat {avg_r:.4g} "
                    f"(good chain mixing during burn-in)"
                )
            elif avg_r < 1.2:
                evidence.append(
                    f"Log file: final avg R-hat {avg_r:.4g} "
                    f"(acceptable chain mixing during burn-in)"
                )
            else:
                evidence.append(
                    f"Log file: final avg R-hat {avg_r:.4g} "
                    f"(poor chain mixing during burn-in)"
                )
                warnings.append(
                    f"High burn-in R-hat ({avg_r:.4g}) suggests poor "
                    f"mixing between chains"
                )

    return _build_result(
        stability_status=status,
        assessment_source=assessment_source,
        iteration_status=iter_status,
        evidence=evidence,
        warnings=warnings,
        details=details,
    )


def _build_result(stability_status, assessment_source, iteration_status,
                  evidence, warnings, details):
    """Construct the assessment result dict.

    Provides both new canonical keys and legacy aliases (``status``,
    ``confidence``) so that callers migrating incrementally continue to work.
    """
    return {
        # ── Canonical keys ──
        "stability_status": stability_status,
        "assessment_source": assessment_source,
        "iteration_status": iteration_status,
        "formal_convergence_status": FormalConvergenceStatus.NOT_FORMALLY_ASSESSED,
        "formal_convergence_reason": _FORMAL_CONVERGENCE_REASON,
        "stability_method": "PIGEAN_MAX_FRACTIONAL_SEM",
        "evidence": evidence,
        "warnings": warnings,
        "details": details,
        # ── Legacy aliases (backward compat) ──
        "status": stability_status,
        "confidence": assessment_source,
    }


# ── JSON writer ───────────────────────────────────────────────────────────

def write_convergence_json(result, output_dir):
    """Write stability assessment to ``convergence.json`` in *output_dir*.

    Converts enums to their string values for JSON serialisation.
    """
    details = result["details"]
    serialisable = {
        "stability_status": result["stability_status"].value,
        "stability_method": result.get("stability_method",
                                       "PIGEAN_MAX_FRACTIONAL_SEM"),
        "assessment_source": result["assessment_source"],
        "max_fractional_sem": details.get("final_max_sem_ratio"),
        "threshold": details.get("max_frac_sem_threshold", 0.01),
        "formal_convergence_status": (
            result["formal_convergence_status"].value
        ),
        "formal_convergence_reason": result.get("formal_convergence_reason"),
        "num_chains": details.get("num_chains"),
        "iteration_status": result["iteration_status"].value,
        "num_gibbs_iter": details.get("num_gibbs_iter"),
        "max_num_iter": details.get("max_num_iter"),
        "num_gibbs_restarts": details.get("num_gibbs_restarts"),
        "evidence": result["evidence"],
        "warnings": result["warnings"],
        "details": details,
    }
    path = os.path.join(output_dir, "convergence.json")
    with open(path, "w") as fh:
        json.dump(serialisable, fh, indent=2)
    return path
