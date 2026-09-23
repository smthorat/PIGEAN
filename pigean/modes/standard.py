"""The locked Phase 1-5 Gibbs workflow."""

from pigean.modes.base import (
    AdvancedMode,
    CompatibilityStatus,
    InterpretationCompatibility,
    StabilityApplicability,
)


class StandardMode(AdvancedMode):
    NAME = "standard"
    ENGINE_SUBCOMMAND = "gibbs"
    DESCRIPTION = "Validated Gibbs-based PIGEAN workflow."
    STABILITY_APPLICABILITY = StabilityApplicability.APPLICABLE
    INTERPRETATION_COMPATIBILITY = (
        InterpretationCompatibility.FULLY_COMPATIBLE
    )
    COMPATIBILITY = {
        "positive-controls": CompatibilityStatus.SUPPORTED,
        "gene-bayes-factor": CompatibilityStatus.SUPPORTED,
        "exome": CompatibilityStatus.SUPPORTED,
        "gwas": CompatibilityStatus.SUPPORTED,
        "gene-z-score": CompatibilityStatus.ENGINE_BLOCKED,
        "gene-percentile": CompatibilityStatus.ENGINE_BLOCKED,
    }

    def build_engine_args(self, config):
        return []
