"""Verified non-outer-Gibbs prior calculation for positive-control input."""

from pigean.modes.base import (
    AdvancedMode,
    CompatibilityStatus,
    InterpretationCompatibility,
    StabilityApplicability,
)


class NaivePriorsMode(AdvancedMode):
    NAME = "naive-priors"
    ENGINE_SUBCOMMAND = "naive_priors"
    DESCRIPTION = (
        "Compute mean-centered gene priors directly from estimated gene-set "
        "effects, without the outer gene-prior Gibbs update loop."
    )
    STABILITY_APPLICABILITY = StabilityApplicability.NOT_APPLICABLE
    INTERPRETATION_COMPATIBILITY = (
        InterpretationCompatibility.PARTIALLY_COMPATIBLE
    )
    COMPATIBILITY = {
        "positive-controls": CompatibilityStatus.SUPPORTED,
        "gene-bayes-factor": CompatibilityStatus.NOT_VERIFIED,
        "exome": CompatibilityStatus.NOT_VERIFIED,
        "gwas": CompatibilityStatus.NOT_VERIFIED,
        "gene-z-score": CompatibilityStatus.ENGINE_BLOCKED,
        "gene-percentile": CompatibilityStatus.ENGINE_BLOCKED,
    }

    def build_engine_args(self, config):
        return []

    def validate_config(self, config):
        issues = super().validate_config(config)
        if config.get("enable_convergence_trace"):
            issues.append(
                "--enable-convergence-trace is incompatible with "
                "mode 'naive-priors' because the outer Gibbs pathway is not run."
            )
        return issues
