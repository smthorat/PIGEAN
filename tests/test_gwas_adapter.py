"""Tests for Phase 4 GWAS adapter and genome-build handling."""
import sys, os, unittest, tempfile, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.adapters import (
    ValidationStatus, get_adapter, SUPPORTED_ANALYSIS_TYPES
)
from pigean.adapters.gwas import GwasAdapter
from pigean.references import normalize_genome_build, GENOME_BUILD_ALIASES

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "gwas")


# ─── GwasAdapter Registry ───────────────────────────────────────────────

class TestGwasAdapterRegistry(unittest.TestCase):
    def test_gwas_in_supported_types(self):
        self.assertIn("gwas", SUPPORTED_ANALYSIS_TYPES)

    def test_get_adapter_returns_gwas(self):
        adapter = get_adapter("gwas")
        self.assertIsInstance(adapter, GwasAdapter)

    def test_gwas_analysis_type(self):
        adapter = GwasAdapter()
        self.assertEqual(adapter.ANALYSIS_TYPE, "gwas")

    def test_gwas_engine_flag(self):
        adapter = GwasAdapter()
        self.assertEqual(adapter.ENGINE_FLAG, "--gwas-in")


# ─── GwasAdapter.validate_file ──────────────────────────────────────────

class TestGwasValidateFile(unittest.TestCase):
    def setUp(self):
        self.adapter = GwasAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_valid_file(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_custom_cols_valid(self):
        path = os.path.join(FIXTURES_DIR, "gwas_custom_cols.tsv")
        params = {
            "gwas_chrom_col": "chromosome",
            "gwas_pos_col": "position",
            "gwas_p_col": "pvalue",
        }
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_specified_column_not_in_header(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        params = {"gwas_chrom_col": "NONEXISTENT"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("NONEXISTENT" in i for i in result.issues))

    def test_multiple_missing_columns(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        params = {"gwas_chrom_col": "FOO", "gwas_pos_col": "BAR"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertGreaterEqual(len(result.issues), 2)

    def test_nonexistent_file(self):
        result = self.adapter.validate_file("/does/not/exist.tsv")
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_empty_file(self):
        path = os.path.join(self.tmpdir, "empty.tsv")
        open(path, "w").close()
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_header_only_file(self):
        path = os.path.join(self.tmpdir, "header_only.tsv")
        with open(path, "w") as f:
            f.write("CHR\tBP\tP\tBETA\tSE\n")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("data row" in i.lower() for i in result.issues))

    def test_column_index_in_range(self):
        """Engine allows 1-based integer column indices."""
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        params = {"gwas_chrom_col": "1"}  # Valid: column 1
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_column_index_out_of_range(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        params = {"gwas_chrom_col": "99"}  # Out of range
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_no_column_overrides_passes(self):
        """Without explicit column overrides, engine auto-detects."""
        path = os.path.join(FIXTURES_DIR, "gwas_missing_cols.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_with_blanks_still_valid(self):
        path = os.path.join(FIXTURES_DIR, "gwas_with_blanks.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)


# ─── GwasAdapter.normalize ──────────────────────────────────────────────

class TestGwasNormalize(unittest.TestCase):
    def setUp(self):
        self.adapter = GwasAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_preserves_original(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertTrue(os.path.exists(result.original_path))
        self.assertIn("gwas_original", result.original_path)

    def test_creates_normalized(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertTrue(os.path.exists(result.normalized_path))
        self.assertIn("gwas.tsv", result.normalized_path)

    def test_row_counts(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertEqual(result.original_rows, 10)
        self.assertEqual(result.normalized_rows, 10)
        self.assertEqual(result.removed_rows, 0)

    def test_blank_lines_removed(self):
        path = os.path.join(FIXTURES_DIR, "gwas_with_blanks.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertEqual(result.original_rows, 5)  # 3 data + 2 blank
        self.assertEqual(result.normalized_rows, 3)
        self.assertEqual(result.blanks_removed, 2)

    def test_transformations_note(self):
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertTrue(len(result.transformations) > 0)
        self.assertTrue(any("No wrapper-level" in t for t in result.transformations))

    def test_normalized_content_matches_original(self):
        """Wrapper does NOT transform GWAS data — just copies."""
        path = os.path.join(FIXTURES_DIR, "valid_gwas_hg19.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        with open(path) as f:
            original = f.read()
        with open(result.normalized_path) as f:
            normalized = f.read()
        self.assertEqual(original, normalized)


# ─── GwasAdapter.build_engine_args ──────────────────────────────────────

class TestGwasBuildEngineArgs(unittest.TestCase):
    def setUp(self):
        self.adapter = GwasAdapter()

    def test_minimal_args(self):
        args = self.adapter.build_engine_args("/tmp/gwas.tsv")
        self.assertEqual(args, ["--gwas-in", "/tmp/gwas.tsv"])

    def test_column_overrides(self):
        params = {
            "gwas_chrom_col": "CHR",
            "gwas_pos_col": "BP",
            "gwas_p_col": "P",
        }
        args = self.adapter.build_engine_args("/tmp/gwas.tsv", params)
        self.assertIn("--gwas-chrom-col", args)
        idx = args.index("--gwas-chrom-col")
        self.assertEqual(args[idx + 1], "CHR")
        self.assertIn("--gwas-pos-col", args)
        self.assertIn("--gwas-p-col", args)

    def test_all_column_overrides(self):
        params = {
            "gwas_chrom_col": "CHR",
            "gwas_pos_col": "BP",
            "gwas_p_col": "P",
            "gwas_beta_col": "BETA",
            "gwas_se_col": "SE",
            "gwas_n_col": "N",
            "gwas_freq_col": "FRQ",
            "gwas_locus_col": "LOC",
            "gwas_filter_col": "FLAG",
            "gwas_filter_value": "PASS",
        }
        args = self.adapter.build_engine_args("/tmp/gwas.tsv", params)
        for flag in ["--gwas-chrom-col", "--gwas-pos-col", "--gwas-p-col",
                      "--gwas-beta-col", "--gwas-se-col", "--gwas-n-col",
                      "--gwas-freq-col", "--gwas-locus-col",
                      "--gwas-filter-col", "--gwas-filter-value"]:
            self.assertIn(flag, args, f"Missing flag: {flag}")

    def test_scalar_gwas_n(self):
        params = {"gwas_n": 50000.0}
        args = self.adapter.build_engine_args("/tmp/gwas.tsv", params)
        self.assertIn("--gwas-n", args)
        idx = args.index("--gwas-n")
        self.assertEqual(args[idx + 1], "50000.0")

    def test_does_not_include_gene_loc_flags(self):
        """Gene location files come from engine.py, not the adapter."""
        params = {"gwas_chrom_col": "CHR"}
        args = self.adapter.build_engine_args("/tmp/gwas.tsv", params)
        self.assertNotIn("--gene-loc-file", args)
        self.assertNotIn("--gene-loc-file-huge", args)
        self.assertNotIn("--exons-loc-file-huge", args)

    def test_empty_params_same_as_none(self):
        args_none = self.adapter.build_engine_args("/tmp/gwas.tsv")
        args_empty = self.adapter.build_engine_args("/tmp/gwas.tsv", {})
        self.assertEqual(args_none, args_empty)


# ─── Genome Build Normalization ─────────────────────────────────────────

class TestGenomeBuildNormalization(unittest.TestCase):
    def test_hg19_passthrough(self):
        self.assertEqual(normalize_genome_build("hg19"), "hg19")

    def test_hg38_passthrough(self):
        self.assertEqual(normalize_genome_build("hg38"), "hg38")

    def test_grch37_to_hg19(self):
        self.assertEqual(normalize_genome_build("GRCh37"), "hg19")

    def test_grch38_to_hg38(self):
        self.assertEqual(normalize_genome_build("GRCh38"), "hg38")

    def test_case_insensitive(self):
        self.assertEqual(normalize_genome_build("HG19"), "hg19")
        self.assertEqual(normalize_genome_build("GRCH37"), "hg19")
        self.assertEqual(normalize_genome_build("HG38"), "hg38")
        self.assertEqual(normalize_genome_build("GRCH38"), "hg38")

    def test_ncbi37_to_hg19(self):
        self.assertEqual(normalize_genome_build("ncbi37"), "hg19")
        self.assertEqual(normalize_genome_build("NCBI37"), "hg19")

    def test_whitespace_stripped(self):
        self.assertEqual(normalize_genome_build(" hg19 "), "hg19")
        self.assertEqual(normalize_genome_build("  GRCh38  "), "hg38")

    def test_unknown_build_raises(self):
        with self.assertRaises(ValueError) as cm:
            normalize_genome_build("hg100")
        self.assertIn("hg100", str(cm.exception))

    def test_empty_string_raises(self):
        with self.assertRaises(ValueError):
            normalize_genome_build("")

    def test_all_aliases_covered(self):
        for alias, expected in GENOME_BUILD_ALIASES.items():
            result = normalize_genome_build(alias)
            self.assertEqual(result, expected,
                             f"Alias '{alias}' should map to '{expected}', got '{result}'")


# ─── Reference Resolution: hg38 files missing ───────────────────────────

class TestHg38ReferencesNotYetAvailable(unittest.TestCase):
    def test_hg38_references_raise_file_not_found(self):
        """hg38 reference files do not exist yet — must fail cleanly."""
        from pigean.references import resolve_reference_paths, resolve_base_dir
        base_dir = resolve_base_dir()
        with self.assertRaises(FileNotFoundError):
            resolve_reference_paths("hg38", base_dir)

    def test_hg19_references_still_resolve(self):
        """hg19 references must still work."""
        from pigean.references import resolve_reference_paths, resolve_base_dir
        base_dir = resolve_base_dir()
        paths = resolve_reference_paths("hg19", base_dir)
        self.assertIn("gene_loc", paths)
        self.assertIn("tss_loc", paths)
        self.assertIn("exons_loc", paths)
        for path in paths.values():
            self.assertTrue(os.path.isfile(path), f"Missing: {path}")


# ─── GWAS in Full Engine Command ────────────────────────────────────────

class TestGwasEngineCommand(unittest.TestCase):
    def setUp(self):
        from pigean.config import DEFAULTS
        self.config = dict(DEFAULTS)
        self.config["analysis"] = "gwas"
        self.config["genome_build"] = "hg19"
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        self.base_dir = repo_root
        self.output_dir = "/tmp/test_gwas_output"
        self.reference_paths = {
            "gene_loc": os.path.join(repo_root, "data", "NCBI37.3.plink.gene.loc"),
            "tss_loc": os.path.join(repo_root, "data", "refGene_hg19_TSS.subset.loc"),
            "exons_loc": os.path.join(repo_root, "data", "NCBI37.3.plink.gene.exons.loc"),
        }
        self.gene_set_paths = [
            os.path.join(repo_root, "data", "gene_set_list_mouse_2024.txt"),
            os.path.join(repo_root, "data", "gene_set_list_msigdb_nohp.txt"),
        ]
        self.gene_map_path = os.path.join(repo_root, "data", "portal_gencode.gene.map")

    def test_gwas_args_in_full_command(self):
        from pigean.engine import build_priors_command
        adapter = GwasAdapter()
        evidence_args = adapter.build_engine_args("/tmp/gwas.tsv")
        cmd = build_priors_command(
            self.config, evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path,
            self.output_dir, self.base_dir,
        )
        self.assertIn("--gwas-in", cmd)
        self.assertNotIn("--positive-controls-in", cmd)

    def test_gwas_command_includes_gene_loc_files(self):
        """Gene location files come from engine.py, not the adapter."""
        from pigean.engine import build_priors_command
        adapter = GwasAdapter()
        evidence_args = adapter.build_engine_args("/tmp/gwas.tsv")
        cmd = build_priors_command(
            self.config, evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path,
            self.output_dir, self.base_dir,
        )
        self.assertIn("--gene-loc-file", cmd)
        self.assertIn("--gene-loc-file-huge", cmd)
        self.assertIn("--exons-loc-file-huge", cmd)

    def test_gwas_with_column_overrides_in_full_command(self):
        from pigean.engine import build_priors_command
        adapter = GwasAdapter()
        params = {
            "gwas_chrom_col": "CHR",
            "gwas_pos_col": "BP",
            "gwas_p_col": "P",
            "gwas_n": 50000,
        }
        evidence_args = adapter.build_engine_args("/tmp/gwas.tsv", params)
        cmd = build_priors_command(
            self.config, evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path,
            self.output_dir, self.base_dir,
        )
        self.assertIn("--gwas-chrom-col", cmd)
        self.assertIn("--gwas-n", cmd)
        # Location files still present
        self.assertIn("--gene-loc-file-huge", cmd)


# ─── Config Integration ─────────────────────────────────────────────────

class TestGwasConfigKeys(unittest.TestCase):
    def test_gwas_keys_in_config_file_keys(self):
        from pigean.config import CONFIG_FILE_KEYS
        gwas_keys = [
            "gwas_chrom_col", "gwas_pos_col", "gwas_p_col",
            "gwas_beta_col", "gwas_se_col", "gwas_n_col",
            "gwas_n", "gwas_freq_col", "gwas_locus_col",
            "gwas_filter_col", "gwas_filter_value",
        ]
        for key in gwas_keys:
            self.assertIn(key, CONFIG_FILE_KEYS, f"Missing from CONFIG_FILE_KEYS: {key}")

    def test_gwas_keys_in_cli_keys(self):
        from pigean.config import CLI_KEYS
        gwas_keys = [
            "gwas_chrom_col", "gwas_pos_col", "gwas_p_col",
            "gwas_beta_col", "gwas_se_col", "gwas_n_col",
            "gwas_n", "gwas_freq_col", "gwas_locus_col",
            "gwas_filter_col", "gwas_filter_value",
        ]
        for key in gwas_keys:
            self.assertIn(key, CLI_KEYS, f"Missing from CLI_KEYS: {key}")

    def test_config_resolution_with_gwas_overrides(self):
        from pigean.config import resolve_config
        args = {
            "analysis": "gwas",
            "genome_build": "hg19",
            "gwas_chrom_col": "CHR",
            "gwas_p_col": "PVALUE",
            "gwas_n": 50000.0,
        }
        config = resolve_config(args)
        self.assertEqual(config["analysis"], "gwas")
        self.assertEqual(config["gwas_chrom_col"], "CHR")
        self.assertEqual(config["gwas_p_col"], "PVALUE")
        self.assertEqual(config["gwas_n"], 50000.0)


if __name__ == "__main__":
    unittest.main()
