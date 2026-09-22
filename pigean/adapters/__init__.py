"""Base types and adapter registry for PIGEAN evidence adapters."""

from enum import Enum
from dataclasses import dataclass, field
from typing import List


class ValidationStatus(Enum):
    """Overall validation outcome."""
    PASS = "PASS"
    PASS_WITH_WARNINGS = "PASS_WITH_WARNINGS"
    FAIL = "FAIL"


@dataclass
class FileValidationResult:
    """Result of file-level input validation."""
    status: ValidationStatus
    issues: List[str] = field(default_factory=list)


@dataclass
class NormalizationResult:
    """Result of input normalization (positive controls / background)."""
    normalized_path: str
    original_count: int = 0
    normalized_count: int = 0
    duplicates_removed: int = 0
    blanks_removed: int = 0
    warnings: List[str] = field(default_factory=list)


@dataclass
class GeneQCResult:
    """Result of gene-level quality control."""
    status: ValidationStatus
    issues: List[str] = field(default_factory=list)
    input_count: int = 0
    recognized_count: int = 0
    unresolved_count: int = 0
    with_coordinates_count: int = 0
    with_annotations_count: int = 0
    per_gene: List[dict] = field(default_factory=list)


@dataclass
class GeneSetPreparationResult:
    """Result of custom gene-set validation + conversion."""
    status: ValidationStatus
    original_paths: List[str] = field(default_factory=list)
    normalized_paths: List[str] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)
    gene_sets_converted: int = 0
    gene_sets_copied: int = 0


@dataclass
class BackgroundQCResult:
    """Result of background gene consistency check."""
    status: ValidationStatus
    issues: List[str] = field(default_factory=list)
    total_count: int = 0
    mapped_count: int = 0
    unmapped_count: int = 0
    input_genes_missing_from_background: List[str] = field(default_factory=list)


class EvidenceAdapter:
    """Legacy base class — kept for backward compatibility with Phase 1/2 code.

    New adapters should extend BaseEvidenceAdapter from pigean.adapters.base.
    """

    def validate_file(self, input_path, params):
        raise NotImplementedError

    def normalize(self, input_path, output_dir, params):
        raise NotImplementedError

    def build_engine_args(self, normalized_path, params):
        raise NotImplementedError


def get_adapter(analysis_type):
    """Return the adapter instance for a given --analysis value.

    Raises ValueError for unknown analysis types.
    """
    from pigean.adapters.positive_controls import PositiveControlsAdapter
    from pigean.adapters.bayes_factor import BayesFactorAdapter
    from pigean.adapters.zscore import ZScoreAdapter
    from pigean.adapters.percentile import PercentileAdapter
    from pigean.adapters.exome import ExomeAdapter
    from pigean.adapters.gwas import GwasAdapter

    registry = {
        "positive-controls": PositiveControlsAdapter,
        "gene-bayes-factor": BayesFactorAdapter,
        "gene-z-score": ZScoreAdapter,
        "gene-percentile": PercentileAdapter,
        "exome": ExomeAdapter,
        "gwas": GwasAdapter,
    }

    if analysis_type not in registry:
        supported = ", ".join(sorted(registry.keys()))
        raise ValueError(
            f"Unknown analysis type: '{analysis_type}'. "
            f"Supported types: {supported}"
        )

    return registry[analysis_type]()


SUPPORTED_ANALYSIS_TYPES = [
    "positive-controls",
    "gene-bayes-factor",
    "gene-z-score",
    "gene-percentile",
    "exome",
    "gwas",
]
