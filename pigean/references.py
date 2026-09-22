"""Reference file path resolution for the PIGEAN wrapper.

Handles both Docker (/app/) and local development environments.
Supports hg19 genome build with explicit hg38 placeholder paths.
Provides genome-build normalization for alias handling.
"""

import os

# Genome build aliases → canonical internal form.
# The wrapper accepts any of these (case-insensitive) and normalizes.
GENOME_BUILD_ALIASES = {
    "hg19": "hg19",
    "grch37": "hg19",
    "ncbi37": "hg19",
    "hg38": "hg38",
    "grch38": "hg38",
}

# Reference file paths relative to the base directory.
# hg38 paths are defined for structural readiness but files do NOT exist yet.
REFERENCE_PATHS = {
    "hg19": {
        "gene_loc": os.path.join("data", "NCBI37.3.plink.gene.loc"),
        "tss_loc": os.path.join("data", "refGene_hg19_TSS.subset.loc"),
        "exons_loc": os.path.join("data", "NCBI37.3.plink.gene.exons.loc"),
    },
    "hg38": {
        "gene_loc": os.path.join("data", "hg38.gene.loc"),
        "tss_loc": os.path.join("data", "hg38.tss.loc"),
        "exons_loc": os.path.join("data", "hg38.exons.loc"),
    },
}

# Annotation gene-set files relative to base directory.
# Order matters: mouse first, msigdb second (matches Phase 0 command).
ANNOTATION_PATHS = {
    "default": [
        os.path.join("data", "gene_set_list_mouse_2024.txt"),
        os.path.join("data", "gene_set_list_msigdb_nohp.txt"),
    ],
    "mouse-only": [
        os.path.join("data", "gene_set_list_mouse_2024.txt"),
    ],
    "msigdb-only": [
        os.path.join("data", "gene_set_list_msigdb_nohp.txt"),
    ],
}

GENE_MAP_PATH = os.path.join("data", "portal_gencode.gene.map")


def normalize_genome_build(build_str):
    """Normalize a genome build string to its canonical form.

    Accepts common aliases (case-insensitive):
        hg19, HG19, GRCh37, grch37, NCBI37 → "hg19"
        hg38, HG38, GRCh38, grch38         → "hg38"

    Parameters
    ----------
    build_str : str
        User-provided genome build identifier.

    Returns
    -------
    str
        Canonical genome build ("hg19" or "hg38").

    Raises
    ------
    ValueError
        If the build string is not recognized.
    """
    canonical = GENOME_BUILD_ALIASES.get(build_str.lower().strip())
    if canonical is None:
        supported = ", ".join(sorted(set(GENOME_BUILD_ALIASES.values())))
        aliases = ", ".join(sorted(GENOME_BUILD_ALIASES.keys()))
        raise ValueError(
            f"Unrecognized genome build: '{build_str}'. "
            f"Supported builds: {supported}. "
            f"Accepted aliases: {aliases}."
        )
    return canonical


def resolve_base_dir():
    """Return the application base directory.

    In Docker: /app (where engine/priors.py lives).
    In local development: the repository root (parent of the pigean/ package).
    """
    if os.path.isfile("/app/engine/priors.py"):
        return "/app"
    # Local dev: pigean/ is in repo root, so go up from this file's directory.
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _verify_file(path, description):
    """Raise FileNotFoundError if path does not exist."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Required {description} not found: {path}"
        )


def resolve_reference_paths(genome_build, base_dir):
    """Resolve genome reference file paths.

    Parameters
    ----------
    genome_build : str
        Genome build identifier ("hg19" or "hg38"). Should be normalized
        via normalize_genome_build() before calling.
    base_dir : str
        Application base directory from resolve_base_dir().

    Returns
    -------
    dict
        Keys: gene_loc, tss_loc, exons_loc. Values: absolute paths.

    Raises
    ------
    ValueError
        If genome_build is not supported.
    FileNotFoundError
        If any required reference file is missing (e.g., hg38 not available).
    """
    if genome_build not in REFERENCE_PATHS:
        raise ValueError(
            f"Unsupported genome build: {genome_build}. "
            f"Supported builds: {list(REFERENCE_PATHS.keys())}"
        )

    refs = {}
    for key, rel_path in REFERENCE_PATHS[genome_build].items():
        abs_path = os.path.join(base_dir, rel_path)
        _verify_file(abs_path, f"{genome_build} {key} file")
        refs[key] = abs_path

    return refs


def resolve_gene_set_paths(gene_set_profile, base_dir,
                            custom_normalized_paths=None,
                            custom_action=None):
    """Resolve annotation gene-set file paths, incorporating custom gene sets.

    Parameters
    ----------
    gene_set_profile : str
        Gene-set profile name ("default", "mouse-only", "msigdb-only", "custom").
    base_dir : str
        Application base directory.
    custom_normalized_paths : list of str or None
        Already-validated-and-converted custom gene-set file paths from
        gene_sets.prepare_custom_gene_sets(). Default None.
    custom_action : str or None
        "replace" or "supplement". Only meaningful when custom_normalized_paths
        is non-empty. Default None.

    Returns
    -------
    list of str
        Ordered list of absolute paths to gene-set files for --X-in flags.

    Raises
    ------
    ValueError
        If gene_set_profile is not supported (and not "custom").
    FileNotFoundError
        If any built-in annotation file is missing.
    """
    custom_paths = custom_normalized_paths or []

    # "custom" profile: use default as base for supplement, custom-only for replace
    if gene_set_profile == "custom":
        if custom_action == "replace" or not custom_paths:
            return list(custom_paths)
        # supplement: default base + custom
        base_paths = _resolve_builtin_paths("default", base_dir)
        return base_paths + list(custom_paths)

    # Built-in profile
    if gene_set_profile not in ANNOTATION_PATHS:
        raise ValueError(
            f"Unsupported gene-set profile: {gene_set_profile}. "
            f"Supported: {list(ANNOTATION_PATHS.keys()) + ['custom']}"
        )

    base_paths = _resolve_builtin_paths(gene_set_profile, base_dir)

    # No custom paths — return base only (Phase 1 behavior)
    if not custom_paths:
        return base_paths

    # With custom paths: replace or supplement
    if custom_action == "replace":
        return list(custom_paths)
    else:
        # supplement (default behavior when custom paths exist)
        return base_paths + list(custom_paths)


def _resolve_builtin_paths(profile_name, base_dir):
    """Resolve and verify built-in annotation paths for a profile."""
    paths = []
    for rel_path in ANNOTATION_PATHS[profile_name]:
        abs_path = os.path.join(base_dir, rel_path)
        _verify_file(abs_path, f"annotation file ({profile_name})")
        paths.append(abs_path)
    return paths


def resolve_gene_map_path(base_dir):
    """Resolve the gene map file path.

    Parameters
    ----------
    base_dir : str
        Application base directory.

    Returns
    -------
    str
        Absolute path to portal_gencode.gene.map.

    Raises
    ------
    FileNotFoundError
        If the gene map file is missing.
    """
    abs_path = os.path.join(base_dir, GENE_MAP_PATH)
    _verify_file(abs_path, "gene map file")
    return abs_path
