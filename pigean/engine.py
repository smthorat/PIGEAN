"""Engine command builder and runner for the PIGEAN wrapper.

Builds the exact priors.py command from resolved configuration and
reference paths, then executes it with full output capture.
"""

import os
import subprocess


def build_priors_command(config, evidence_engine_args, reference_paths,
                        gene_set_paths, gene_map_path, output_dir, base_dir,
                        background_normalized_path=None,
                        enable_convergence_trace=False,
                        engine_subcommand="gibbs",
                        mode_engine_args=None):
    """Build the priors.py command as a Python list.

    Uses full engine option names (not optparse prefix abbreviations).

    Parameters
    ----------
    config : dict
        Resolved configuration from config.resolve_config().
    evidence_engine_args : list of str
        Evidence-specific engine arguments from the adapter's
        build_engine_args() method. For positive-controls:
        ["--positive-controls-in", path]. For gene-bfs:
        ["--gene-bfs-in", path]. Etc.
    reference_paths : dict
        Genome reference paths from references.resolve_reference_paths().
    gene_set_paths : list of str
        Ordered annotation file paths from references.resolve_gene_set_paths().
    gene_map_path : str
        Path to the gene map from references.resolve_gene_map_path().
    output_dir : str
        Output directory for engine results.
    base_dir : str
        Application base directory (contains engine/priors.py).
    background_normalized_path : str or None
        Path to the normalized background gene list. When provided,
        adds --positive-controls-all-in to the engine command.

    Returns
    -------
    list of str
        The command as a list suitable for subprocess.
    """
    cmd = [
        "python3", "-u", os.path.join(base_dir, "engine", "priors.py"),
        engine_subcommand,
    ]

    # Mode-specific options are deliberately small and independently
    # validated by pigean.modes. Evidence options are appended separately.
    cmd.extend(mode_engine_args or [])

    # Annotation gene-set files (order matters — matches Phase 0)
    for gs_path in gene_set_paths:
        cmd.extend(["--X-in", gs_path])

    # Gene map
    cmd.extend(["--gene-map-in", gene_map_path])

    # Max gene sets
    cmd.extend(["--max-num-gene-sets", str(config["max_num_gene_sets"])])

    # Output files (using full flag names)
    cmd.extend(["--gene-stats-out", os.path.join(output_dir, "gs.out")])
    cmd.extend(["--gene-set-stats-out", os.path.join(output_dir, "gss.out")])
    cmd.extend(["--gene-gene-set-stats-out",
                 os.path.join(output_dir, "ggss.out")])
    cmd.extend(["--params-out", os.path.join(output_dir, "p.out")])

    # Debug level
    cmd.extend(["--debug-level", str(config["debug_level"])])

    # Evidence-specific arguments (from adapter)
    cmd.extend(evidence_engine_args)

    # Genome reference files
    cmd.extend(["--gene-loc-file", reference_paths["gene_loc"]])
    cmd.extend(["--gene-loc-file-huge", reference_paths["tss_loc"]])
    cmd.extend(["--exons-loc-file-huge", reference_paths["exons_loc"]])

    # Filter values
    cmd.extend(["--gene-filter-value", str(config["gene_filter_value"])])
    cmd.extend(["--gene-set-filter-value",
                 str(config["gene_set_filter_value"])])

    # Background genes (eligible gene population)
    if background_normalized_path is not None:
        cmd.extend(["--positive-controls-all-in",
                     background_normalized_path])

    # Phase 5: Optional convergence trace output
    if enable_convergence_trace:
        cmd.extend(["--gene-set-stats-trace-out",
                     os.path.join(output_dir, "gss_trace.out")])

    return cmd


def save_command(cmd, output_dir):
    """Write priors_command.txt — the exact engine invocation."""
    path = os.path.join(output_dir, "priors_command.txt")

    header = " ".join(cmd[:4])
    body_lines = []
    i = 4
    while i < len(cmd):
        token = cmd[i]
        if token.startswith("--") and i + 1 < len(cmd) and not cmd[i + 1].startswith("--"):
            body_lines.append(f"    {token} {cmd[i + 1]}")
            i += 2
        elif token.startswith("--"):
            body_lines.append(f"    {token}")
            i += 1
        else:
            body_lines.append(f"    {token}")
            i += 1

    with open(path, "w") as f:
        if body_lines:
            f.write(header + " \\\n")
            for j, line in enumerate(body_lines):
                if j < len(body_lines) - 1:
                    f.write(line + " \\\n")
                else:
                    f.write(line + "\n")
        else:
            f.write(header + "\n")

    return path


def run_engine(cmd, output_dir):
    """Execute the priors.py engine command.

    Captures stdout and stderr to pigean_run.log in the output directory.
    Returns the process exit code.
    """
    log_path = os.path.join(output_dir, "pigean_run.log")

    with open(log_path, "w") as log_file:
        result = subprocess.run(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )

    return result.returncode
