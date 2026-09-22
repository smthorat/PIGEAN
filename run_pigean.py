#!/usr/bin/env python3
"""PIGEAN Pipeline Wrapper — user-facing interface to priors.py."""

import argparse
import os
import sys
import shutil
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pigean import __version__
from pigean.config import resolve_config, write_resolved_config
from pigean.references import (
    resolve_base_dir, resolve_reference_paths,
    resolve_gene_set_paths, resolve_gene_map_path,
    normalize_genome_build,
)
from pigean.adapters import (
    ValidationStatus, SUPPORTED_ANALYSIS_TYPES, get_adapter
)
from pigean.adapters.positive_controls import PositiveControlsAdapter
from pigean.adapters.gene_sets import prepare_custom_gene_sets
from pigean.validation import (
    validate_input_genes, write_gene_qc_tsv,
    compute_validation_status, validate_background,
    validate_evidence_genes, write_evidence_qc_tsv,
)
from pigean.engine import build_priors_command, save_command, run_engine
from pigean.manifest import generate_manifest, write_manifest
from pigean.report import generate_report, generate_report_html
from pigean.parsers import (
    parse_gene_stats, parse_gene_set_stats,
    parse_gene_gene_set_stats, parse_params,
)
from pigean.convergence import assess_convergence, write_convergence_json
from pigean.interpretation import interpret_results


def _write_fail_manifest(config, base_dir, output_dir, start_time,
                         file_result=None, norm_result=None, qc_result=None,
                         validation_issues=None, **extra_kwargs):
    """Write a manifest for a failed validation and exit."""
    end_time = datetime.now()
    manifest = generate_manifest(
        config=config,
        file_result=file_result,
        norm_result=norm_result,
        qc_result=qc_result,
        validation_status=ValidationStatus.FAIL,
        validation_issues=validation_issues or [],
        execution_status="NOT_RUN",
        exit_code=None,
        base_dir=base_dir,
        reference_paths=None,
        gene_set_paths=None,
        gene_map_path=None,
        output_dir=output_dir,
        start_time=start_time,
        end_time=end_time,
        **extra_kwargs,
    )
    write_manifest(manifest, output_dir)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=f"PIGEAN Pipeline Wrapper v{__version__}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Positive controls (default)
  python3 run_pigean.py --input ex/gene_list --output runs/my_run

  # Gene-level Bayes factors
  python3 run_pigean.py --analysis gene-bayes-factor \\
      --input evidence.tsv --output runs/bf_run

  # Gene-level scores (Z-score/log-odds mode)
  python3 run_pigean.py --analysis gene-z-score \\
      --input scores.tsv --gene-column Gene --score-column Score \\
      --output runs/zscore_run

  # Exome associations
  python3 run_pigean.py --analysis exome \\
      --input exome_results.tsv --output runs/exome_run
        """,
    )

    # Required
    parser.add_argument("--input", required=True,
                        help="Input evidence file")
    parser.add_argument("--output", required=True,
                        help="Output directory for results")

    # Analysis settings
    parser.add_argument("--analysis", default=None,
                        choices=SUPPORTED_ANALYSIS_TYPES,
                        help="Analysis/evidence type (default: positive-controls)")
    parser.add_argument("--gene-sets", default=None,
                        choices=["default", "mouse-only", "msigdb-only", "custom"],
                        help="Gene set profile (default: default)")
    parser.add_argument("--genome-build", default=None,
                        help="Genome build (default: hg19). "
                             "REQUIRED for GWAS mode. "
                             "Accepts: hg19, GRCh37, hg38, GRCh38.")
    parser.add_argument("--preset", default=None,
                        choices=["standard"],
                        help="Parameter preset (default: standard)")

    # Optional
    parser.add_argument("--config", default=None,
                        help="Optional JSON config file for parameter overrides")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing output directory")

    # Phase 2: Custom gene sets
    parser.add_argument("--custom-gene-set-files", nargs="+", default=None,
                        help="One or more custom gene-set files")
    parser.add_argument("--custom-gene-set-format", default=None,
                        choices=["engine", "gmt"],
                        help="Format of custom gene-set files")
    parser.add_argument("--custom-gene-set-action", default=None,
                        choices=["replace", "supplement"],
                        help="Whether custom gene sets replace or supplement the base profile")

    # Phase 2: Background genes
    parser.add_argument("--background", default=None,
                        help="Background gene list file for --positive-controls-all-in")

    # Phase 3: Column mapping (for tabular evidence files)
    parser.add_argument("--gene-column", default=None,
                        help="Column name for gene identifiers in evidence file")
    parser.add_argument("--score-column", default=None,
                        help="Column name for evidence scores in evidence file")

    # Phase 3: Percentile-specific
    parser.add_argument("--higher-is-better", action="store_true", default=None,
                        help="For gene-percentile: higher score = stronger evidence")

    # Phase 3: Exome-specific column overrides
    parser.add_argument("--exomes-gene-col", default=None,
                        help="Gene column name for exome input")
    parser.add_argument("--exomes-p-col", default=None,
                        help="P-value column name for exome input")
    parser.add_argument("--exomes-beta-col", default=None,
                        help="Beta/effect column name for exome input")
    parser.add_argument("--exomes-se-col", default=None,
                        help="Standard error column name for exome input")
    parser.add_argument("--exomes-n-col", default=None,
                        help="Sample size column name for exome input")
    parser.add_argument("--exomes-n", type=float, default=None,
                        help="Global sample size for exome input")

    # Phase 4: GWAS column overrides
    parser.add_argument("--gwas-chrom-col", default=None,
                        help="Chromosome column name in GWAS file")
    parser.add_argument("--gwas-pos-col", default=None,
                        help="Position column name in GWAS file")
    parser.add_argument("--gwas-p-col", default=None,
                        help="P-value column name in GWAS file")
    parser.add_argument("--gwas-beta-col", default=None,
                        help="Beta/effect size column name in GWAS file")
    parser.add_argument("--gwas-se-col", default=None,
                        help="Standard error column name in GWAS file")
    parser.add_argument("--gwas-n-col", default=None,
                        help="Sample size column name in GWAS file")
    parser.add_argument("--gwas-n", type=float, default=None,
                        help="Global sample size for GWAS (scalar)")
    parser.add_argument("--gwas-freq-col", default=None,
                        help="Allele frequency column name in GWAS file")
    parser.add_argument("--gwas-locus-col", default=None,
                        help="Locus column (chr:pos compound) in GWAS file")
    parser.add_argument("--gwas-filter-col", default=None,
                        help="Row-filter column name in GWAS file")
    parser.add_argument("--gwas-filter-value", default=None,
                        help="Row-filter match value for GWAS file")

    # Phase 5: Convergence diagnostics
    parser.add_argument("--enable-convergence-trace", action="store_true",
                        default=None,
                        help="Enable convergence trace output (gss_trace.out)")

    args = parser.parse_args()

    start_time = datetime.now()

    # ── Step 1: Create clean output directory ──
    output_dir = os.path.abspath(args.output)

    if os.path.exists(output_dir):
        if os.listdir(output_dir):
            if not args.overwrite:
                print(f"ERROR: Output directory is not empty: {output_dir}")
                print("Use --overwrite to replace existing contents.")
                sys.exit(1)
            shutil.rmtree(output_dir)

    os.makedirs(output_dir, exist_ok=True)

    # ── Step 2: Resolve configuration ──
    config = resolve_config(args, config_file=args.config)

    # ── Step 2a: Normalize and validate genome build ──
    analysis_type = config["analysis"]
    is_gwas = (analysis_type == "gwas")

    # Normalize genome build aliases (hg19/GRCh37/etc.)
    if config.get("genome_build"):
        try:
            config["genome_build"] = normalize_genome_build(config["genome_build"])
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # GWAS mode requires explicit genome build — never silently default
    if is_gwas and args.genome_build is None:
        print("ERROR: GWAS analysis requires --genome-build to be specified explicitly.")
        print("       The genome build determines which reference files are used for")
        print("       SNP-to-gene mapping. Mixing builds silently produces wrong results.")
        print("       Example: --genome-build hg19  (or GRCh37, hg38, GRCh38)")
        sys.exit(1)

    # ── Step 2b: Resolve base directory ──
    base_dir = resolve_base_dir()

    # ── Step 3: Write resolved config ──
    write_resolved_config(config, output_dir)

    # ── Step 4: Select evidence adapter ──
    try:
        adapter = get_adapter(analysis_type)
    except ValueError as e:
        print(f"ERROR: {e}")
        _write_fail_manifest(config, base_dir, output_dir, start_time,
                             validation_issues=[str(e)])

    is_positive_controls = (analysis_type == "positive-controls")

    # ── Step 5: Validate primary input file ──
    input_path = os.path.abspath(args.input)
    file_result = adapter.validate_file(input_path, config)

    if file_result.status == ValidationStatus.FAIL:
        print(f"ERROR: Input validation failed:")
        for issue in file_result.issues:
            print(f"  - {issue}")
        _write_fail_manifest(config, base_dir, output_dir, start_time,
                             file_result=file_result,
                             validation_issues=file_result.issues)

    # ── Step 6: Preserve original + normalize ──
    if is_positive_controls:
        norm_result = adapter.normalize(input_path, output_dir, config)
        normalized_path = norm_result.normalized_path
        evidence_norm_result = None
    else:
        norm_result = None
        evidence_norm_result = adapter.normalize(input_path, output_dir, config)
        normalized_path = evidence_norm_result.normalized_path

    # ── Step 6b: Background file (if requested) ──
    background_norm_result = None
    background_file_result = None
    background_qc_result = None
    if args.background:
        bg_path = os.path.abspath(args.background)
        bg_adapter = PositiveControlsAdapter()
        background_file_result = bg_adapter.validate_file(bg_path, config)

        if background_file_result.status == ValidationStatus.FAIL:
            print(f"ERROR: Background file validation failed:")
            for issue in background_file_result.issues:
                print(f"  - {issue}")
            _write_fail_manifest(config, base_dir, output_dir, start_time,
                                 file_result=file_result, norm_result=norm_result,
                                 validation_issues=background_file_result.issues)

        background_norm_result = bg_adapter.normalize(
            bg_path, output_dir, config,
            input_label="background_original",
            normalized_label="background",
        )

    # ── Step 6c: Validate + prepare custom gene sets ──
    custom_gs_prep = prepare_custom_gene_sets(
        custom_files=[os.path.abspath(f) for f in (args.custom_gene_set_files or [])],
        custom_format=config.get("custom_gene_set_format"),
        custom_action=config.get("custom_gene_set_action"),
        gene_set_profile=config["gene_sets"],
        output_dir=output_dir,
    )
    if custom_gs_prep.status == ValidationStatus.FAIL:
        print(f"ERROR: Custom gene-set preparation failed:")
        for issue in custom_gs_prep.issues:
            print(f"  - {issue}")
        _write_fail_manifest(config, base_dir, output_dir, start_time,
                             file_result=file_result, norm_result=norm_result,
                             validation_issues=custom_gs_prep.issues)

    # ── Step 7: Resolve references ──
    try:
        reference_paths = resolve_reference_paths(config["genome_build"], base_dir)
        gene_set_paths = resolve_gene_set_paths(
            config["gene_sets"], base_dir,
            custom_normalized_paths=custom_gs_prep.normalized_paths,
            custom_action=config.get("custom_gene_set_action"),
        )
        gene_map_path = resolve_gene_map_path(base_dir)
    except FileNotFoundError as e:
        print(f"ERROR: Required reference file missing: {e}")
        _write_fail_manifest(config, base_dir, output_dir, start_time,
                             file_result=file_result, norm_result=norm_result,
                             validation_issues=[f"Missing reference: {e}"])

    # ── Step 8: Gene/reference QC ──
    if is_positive_controls:
        qc_result = validate_input_genes(
            normalized_path, gene_map_path,
            reference_paths["gene_loc"], gene_set_paths
        )
        write_gene_qc_tsv(qc_result, output_dir)
        evidence_qc_result = None
    elif is_gwas:
        # GWAS is SNP-level — gene QC does not apply.
        # The engine maps SNPs to genes internally using reference files.
        # Validation was done at file level in Step 5.
        from pigean.adapters.base import EvidenceQCResult
        evidence_qc_result = EvidenceQCResult(
            evidence_type="gwas",
            status=ValidationStatus.PASS,
            issues=[],
            total_genes=0,
            recognized_genes=0,
        )
        evidence_qc_result.issues.append(
            "GWAS input is SNP-level; gene QC not applicable. "
            "Engine performs SNP-to-gene mapping using reference files."
        )
        evidence_qc_result.status = ValidationStatus.PASS_WITH_WARNINGS
        write_evidence_qc_tsv(evidence_qc_result, output_dir)
        qc_result = evidence_qc_result
    else:
        evidence_qc_result = validate_evidence_genes(
            normalized_path, gene_map_path, config
        )
        write_evidence_qc_tsv(evidence_qc_result, output_dir)
        qc_result = evidence_qc_result

    # ── Step 8b: Background QC ──
    if background_norm_result is not None:
        background_qc_result = validate_background(
            background_norm_result.normalized_path,
            normalized_path if is_positive_controls else None,
            gene_map_path,
        )

    # ── Step 9: Combined validation status ──
    final_status = compute_validation_status(
        file_result, qc_result,
        background_file_result=background_file_result,
        background_qc_result=background_qc_result,
        custom_gs_prep_result=custom_gs_prep,
    )
    all_issues = file_result.issues + qc_result.issues
    if background_file_result and background_file_result.issues:
        all_issues += background_file_result.issues
    if background_qc_result and background_qc_result.issues:
        all_issues += background_qc_result.issues
    if custom_gs_prep.issues:
        all_issues += custom_gs_prep.issues

    if final_status == ValidationStatus.PASS_WITH_WARNINGS:
        print("WARNINGS:")
        for issue in all_issues:
            print(f"  - {issue}")

    # ── Step 10: If FAIL, write manifest and exit ──
    if final_status == ValidationStatus.FAIL:
        print(f"ERROR: Validation failed — engine will NOT run:")
        for issue in all_issues:
            print(f"  - {issue}")

        end_time = datetime.now()
        manifest = generate_manifest(
            config=config,
            file_result=file_result,
            norm_result=norm_result,
            qc_result=qc_result,
            validation_status=final_status,
            validation_issues=all_issues,
            execution_status="NOT_RUN",
            exit_code=None,
            base_dir=base_dir,
            reference_paths=reference_paths,
            gene_set_paths=gene_set_paths,
            gene_map_path=gene_map_path,
            output_dir=output_dir,
            start_time=start_time,
            end_time=end_time,
            evidence_norm_result=evidence_norm_result,
        )
        write_manifest(manifest, output_dir)
        sys.exit(1)

    # ── Steps 11-12: Build and save priors.py command ──
    evidence_engine_args = adapter.build_engine_args(normalized_path, config)

    cmd = build_priors_command(
        config, evidence_engine_args, reference_paths,
        gene_set_paths, gene_map_path, output_dir, base_dir,
        background_normalized_path=(
            background_norm_result.normalized_path
            if background_norm_result else None
        ),
        enable_convergence_trace=bool(config.get("enable_convergence_trace")),
    )
    save_command(cmd, output_dir)

    # ── Steps 13-14: Run priors.py ──
    print(f"Running PIGEAN engine ({analysis_type})...")
    exit_code = run_engine(cmd, output_dir)

    execution_status = "COMPLETED" if exit_code == 0 else "FAILED"
    if exit_code != 0:
        print(f"WARNING: Engine exited with code {exit_code}")
    else:
        print(f"Engine completed successfully.")

    # ── Phase 5: Parse outputs, assess convergence, interpret ──
    parsed_outputs = None
    convergence_result = None
    interpretation_result = None

    if exit_code == 0:
        parsed_gs = parse_gene_stats(os.path.join(output_dir, "gs.out"))
        parsed_gss = parse_gene_set_stats(os.path.join(output_dir, "gss.out"))
        parsed_ggss = parse_gene_gene_set_stats(
            os.path.join(output_dir, "ggss.out"))
        parsed_params = parse_params(os.path.join(output_dir, "p.out"))

        parsed_outputs = {
            "gs": parsed_gs, "gss": parsed_gss,
            "ggss": parsed_ggss, "p": parsed_params,
        }

        convergence_result = assess_convergence(
            parsed_params,
            log_path=os.path.join(output_dir, "pigean_run.log"),
        )
        write_convergence_json(convergence_result, output_dir)
        print(f"PIGEAN stability: "
              f"{convergence_result['stability_status'].value} "
              f"({convergence_result['assessment_source']})")

        interpretation_result = interpret_results(
            parsed_gs, parsed_gss, parsed_ggss, parsed_params,
            convergence_result, analysis_type,
        )
    else:
        # Engine failed — assess what we can
        parsed_params = parse_params(os.path.join(output_dir, "p.out"))
        if parsed_params is not None:
            convergence_result = assess_convergence(
                parsed_params,
                log_path=os.path.join(output_dir, "pigean_run.log"),
            )
            write_convergence_json(convergence_result, output_dir)

    # ── Step 15: Write manifest ──
    end_time = datetime.now()
    manifest = generate_manifest(
        config=config,
        file_result=file_result,
        norm_result=norm_result,
        qc_result=qc_result,
        validation_status=final_status,
        validation_issues=all_issues,
        execution_status=execution_status,
        exit_code=exit_code,
        base_dir=base_dir,
        reference_paths=reference_paths,
        gene_set_paths=gene_set_paths,
        gene_map_path=gene_map_path,
        output_dir=output_dir,
        start_time=start_time,
        end_time=end_time,
        background_path=(os.path.abspath(args.background)
                         if args.background else None),
        background_norm_result=background_norm_result,
        background_qc_result=background_qc_result,
        custom_gs_prep_result=custom_gs_prep,
        evidence_norm_result=evidence_norm_result,
        parsed_outputs=parsed_outputs,
        convergence_result=convergence_result,
        interpretation_result=interpretation_result,
    )
    write_manifest(manifest, output_dir)

    # ── Step 16: Generate report ──
    if exit_code == 0:
        generate_report(config, qc_result, manifest, output_dir,
                        evidence_norm_result=evidence_norm_result,
                        interpretation_result=interpretation_result,
                        convergence_result=convergence_result)
        generate_report_html(config, manifest, output_dir,
                             interpretation_result=interpretation_result,
                             convergence_result=convergence_result)
        print(f"Report written to {os.path.join(output_dir, 'report.txt')}")
        print(f"HTML report written to {os.path.join(output_dir, 'report.html')}")

    print(f"\nResults in: {output_dir}")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
