"""High-level advanced-mode contract for the PIGEAN wrapper.

Modes select the engine's statistical pathway.  They are intentionally
separate from evidence adapters, which describe the user's input evidence.
"""

from abc import ABC, abstractmethod
from enum import Enum


class CompatibilityStatus(Enum):
    SUPPORTED = "SUPPORTED"
    ENGINE_BLOCKED = "ENGINE_BLOCKED"
    INCOMPATIBLE = "INCOMPATIBLE"
    NOT_VERIFIED = "NOT_VERIFIED"


class InterpretationCompatibility(Enum):
    FULLY_COMPATIBLE = "FULLY_COMPATIBLE"
    PARTIALLY_COMPATIBLE = "PARTIALLY_COMPATIBLE"
    MODE_SPECIFIC_INTERPRETATION_REQUIRED = (
        "MODE_SPECIFIC_INTERPRETATION_REQUIRED"
    )
    NOT_INTERPRETABLE_WITH_CURRENT_CONTRACT = (
        "NOT_INTERPRETABLE_WITH_CURRENT_CONTRACT"
    )


class StabilityApplicability(Enum):
    APPLICABLE = "APPLICABLE"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class AdvancedMode(ABC):
    """Small interface implemented by every researcher-facing model mode."""

    NAME = NotImplemented
    ENGINE_SUBCOMMAND = NotImplemented
    DESCRIPTION = NotImplemented
    STABILITY_APPLICABILITY = StabilityApplicability.APPLICABLE
    INTERPRETATION_COMPATIBILITY = (
        InterpretationCompatibility.FULLY_COMPATIBLE
    )
    COMPATIBILITY = {}
    UNEXPOSED_ADVANCED_PARAMETERS = ("anchor", "phi", "alpha0")

    def compatibility_for(self, evidence_type):
        return self.COMPATIBILITY.get(
            evidence_type, CompatibilityStatus.NOT_VERIFIED
        )

    def compatible_evidence_types(self):
        return [
            evidence
            for evidence, status in self.COMPATIBILITY.items()
            if status == CompatibilityStatus.SUPPORTED
        ]

    def compatible_gene_set_profiles(self):
        return ["default", "mouse-only", "msigdb-only", "custom"]

    def validate_config(self, config):
        """Return validation issues; an empty list means the mode is valid."""
        supplied_internal = [
            name for name in self.UNEXPOSED_ADVANCED_PARAMETERS
            if config.get(name) is not None
        ]
        issues = []
        if supplied_internal:
            issues.append(
                "Internal advanced parameter(s) are not exposed by the "
                "validated wrapper: " + ", ".join(supplied_internal) + "."
            )
        evidence_type = config.get("analysis", "positive-controls")
        status = self.compatibility_for(evidence_type)
        if status == CompatibilityStatus.SUPPORTED:
            return issues
        issues.append(
            f"Mode '{self.NAME}' with evidence type '{evidence_type}' is "
            f"{status.value}. The engine will not be run."
        )
        return issues

    @abstractmethod
    def build_engine_args(self, config):
        """Return mode-specific engine arguments after the subcommand."""
        raise NotImplementedError

    def describe(self, config=None):
        config = config or {}
        evidence_type = config.get("analysis", "positive-controls")
        supplied = {
            name: config[name]
            for name in self.UNEXPOSED_ADVANCED_PARAMETERS
            if config.get(name) is not None
        }
        return {
            "name": self.NAME,
            "description": self.DESCRIPTION,
            "engine_subcommand": self.ENGINE_SUBCOMMAND,
            "user_supplied_parameters": supplied,
            "resolved_parameters": {},
            "engine_arguments": [self.ENGINE_SUBCOMMAND]
            + self.build_engine_args(config),
            "compatibility_status": self.compatibility_for(
                evidence_type
            ).value,
            "stability_applicability": self.STABILITY_APPLICABILITY.value,
            "interpretation_compatibility": (
                self.INTERPRETATION_COMPATIBILITY.value
            ),
        }
