"""Registry for verified high-level PIGEAN model modes."""

from pigean.modes.base import (
    AdvancedMode,
    CompatibilityStatus,
    InterpretationCompatibility,
    StabilityApplicability,
)
from pigean.modes.naive_priors import NaivePriorsMode
from pigean.modes.standard import StandardMode


MODES = {
    "standard": StandardMode,
    "naive-priors": NaivePriorsMode,
}

SUPPORTED_MODES = list(MODES.keys())

COMPATIBILITY_MATRIX = {
    evidence: {
        mode_name: mode_cls().compatibility_for(evidence).value
        for mode_name, mode_cls in MODES.items()
    }
    for evidence in [
        "positive-controls",
        "gene-bayes-factor",
        "exome",
        "gwas",
        "gene-z-score",
        "gene-percentile",
    ]
}


def get_mode(mode_name):
    """Return a verified mode instance; never silently fall back."""
    mode_cls = MODES.get(mode_name)
    if mode_cls is None:
        supported = ", ".join(SUPPORTED_MODES)
        raise ValueError(
            f"Unknown or unexposed mode: '{mode_name}'. "
            f"Supported modes: {supported}."
        )
    return mode_cls()


__all__ = [
    "AdvancedMode",
    "CompatibilityStatus",
    "InterpretationCompatibility",
    "StabilityApplicability",
    "MODES",
    "SUPPORTED_MODES",
    "COMPATIBILITY_MATRIX",
    "get_mode",
]
