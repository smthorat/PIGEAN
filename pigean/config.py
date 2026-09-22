"""Configuration defaults and resolution for the PIGEAN wrapper.

Configuration precedence (later overrides earlier):
    DEFAULTS → preset overrides → JSON config file → explicit CLI args
"""

import json
import os

# Verified against the locked Phase 0 command and engine behavior.
# Each value here matches the effective scientific setting used in Phase 0.
DEFAULTS = {
    "analysis": "positive-controls",
    "preset": "standard",
    "genome_build": "hg19",
    "gene_sets": "default",
    "max_num_gene_sets": 5000,
    "gene_filter_value": 1,
    "gene_set_filter_value": 0.01,
    "debug_level": 3,
}

PRESETS = {
    "standard": {},
}

# Keys that are valid in a user-supplied JSON config file.
CONFIG_FILE_KEYS = {
    "analysis", "preset", "genome_build", "gene_sets",
    "max_num_gene_sets", "gene_filter_value", "gene_set_filter_value",
    "debug_level",
    "custom_gene_set_format", "custom_gene_set_action",
    # Phase 3: evidence-related settings
    "gene_column", "score_column", "higher_is_better",
    "exomes_gene_col", "exomes_p_col", "exomes_beta_col",
    "exomes_se_col", "exomes_n_col", "exomes_n",
    # Phase 4: GWAS column overrides
    "gwas_chrom_col", "gwas_pos_col", "gwas_p_col", "gwas_beta_col",
    "gwas_se_col", "gwas_n_col", "gwas_n", "gwas_freq_col",
    "gwas_locus_col", "gwas_filter_col", "gwas_filter_value",
    # Phase 5: convergence diagnostics
    "enable_convergence_trace",
}

# Keys that correspond to argparse CLI options (with default=None).
CLI_KEYS = {
    "analysis", "preset", "genome_build", "gene_sets",
    "custom_gene_set_format", "custom_gene_set_action",
    # Phase 3
    "gene_column", "score_column", "higher_is_better",
    "exomes_gene_col", "exomes_p_col", "exomes_beta_col",
    "exomes_se_col", "exomes_n_col", "exomes_n",
    # Phase 4: GWAS column overrides
    "gwas_chrom_col", "gwas_pos_col", "gwas_p_col", "gwas_beta_col",
    "gwas_se_col", "gwas_n_col", "gwas_n", "gwas_freq_col",
    "gwas_locus_col", "gwas_filter_col", "gwas_filter_value",
    # Phase 5: convergence diagnostics
    "enable_convergence_trace",
}


def resolve_config(args, config_file=None):
    """Merge configuration from all sources.

    Parameters
    ----------
    args : argparse.Namespace or dict
        Parsed CLI arguments. Only non-None values override.
    config_file : str or None
        Path to a JSON config file (optional).

    Returns
    -------
    dict
        Fully resolved configuration.
    """
    config = dict(DEFAULTS)

    if isinstance(args, dict):
        preset_name = args.get("preset") or config["preset"]
    else:
        preset_name = getattr(args, "preset", None) or config["preset"]

    if preset_name in PRESETS:
        config.update(PRESETS[preset_name])
    config["preset"] = preset_name

    if config_file is not None:
        with open(config_file, "r") as f:
            file_config = json.load(f)
        for key, value in file_config.items():
            if key in CONFIG_FILE_KEYS and value is not None:
                config[key] = value

    if isinstance(args, dict):
        args_dict = args
    else:
        args_dict = vars(args)

    for key in CONFIG_FILE_KEYS:
        value = args_dict.get(key)
        if value is not None:
            config[key] = value

    return config


def write_resolved_config(config, output_dir):
    """Write resolved_config.json to the output directory."""
    path = os.path.join(output_dir, "resolved_config.json")
    with open(path, "w") as f:
        json.dump(config, f, indent=2, sort_keys=True)
    return path
