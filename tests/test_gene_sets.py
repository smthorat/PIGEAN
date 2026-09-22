"""Tests for pigean/adapters/gene_sets.py — GMT conversion, legal combinations,
resolve_gene_set_paths with custom, and background validation."""
import sys, os, tempfile, shutil, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.adapters import (
    ValidationStatus, FileValidationResult, GeneQCResult,
    GeneSetPreparationResult, BackgroundQCResult,
)
from pigean.adapters.gene_sets import (
    validate_gene_set_file, convert_gmt_to_engine_format,
    prepare_custom_gene_sets,
)
from pigean.references import resolve_gene_set_paths

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GENE_MAP = os.path.join(REPO_ROOT, "data", "portal_gencode.gene.map")


# ── GMT Conversion ──

class TestGMTConversion(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_valid_gmt_conversion(self):
        """GMT with name+desc+genes → engine with name+genes."""
        gmt = os.path.join(self.tmpdir, "test.gmt")
        out = os.path.join(self.tmpdir, "test.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATHWAY1\tsome description\tGENE1\tGENE2\tGENE3\n")
        result = convert_gmt_to_engine_format(gmt, out)
        self.assertEqual(result["gene_sets_converted"], 1)
        with open(out) as f:
            line = f.read().strip()
        self.assertEqual(line, "PATHWAY1\tGENE1\tGENE2\tGENE3")

    def test_multi_line_gmt(self):
        """Multiple GMT lines convert correctly."""
        gmt = os.path.join(self.tmpdir, "multi.gmt")
        out = os.path.join(self.tmpdir, "multi.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATH1\tdesc1\tA\tB\n")
            f.write("PATH2\tdesc2\tC\tD\tE\n")
        result = convert_gmt_to_engine_format(gmt, out)
        self.assertEqual(result["gene_sets_converted"], 2)
        with open(out) as f:
            lines = [l.strip() for l in f if l.strip()]
        self.assertEqual(lines[0], "PATH1\tA\tB")
        self.assertEqual(lines[1], "PATH2\tC\tD\tE")

    def test_gmt_missing_genes_raises(self):
        """GMT line with only name+desc (no genes) → ValueError."""
        gmt = os.path.join(self.tmpdir, "bad.gmt")
        out = os.path.join(self.tmpdir, "bad.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATHWAY1\tsome description\n")
        with self.assertRaises(ValueError):
            convert_gmt_to_engine_format(gmt, out)

    def test_gmt_single_column_raises(self):
        """GMT line with only name → ValueError."""
        gmt = os.path.join(self.tmpdir, "single.gmt")
        out = os.path.join(self.tmpdir, "single.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATHWAY1\n")
        with self.assertRaises(ValueError):
            convert_gmt_to_engine_format(gmt, out)

    def test_blank_lines_skipped(self):
        """Blank lines in GMT are ignored."""
        gmt = os.path.join(self.tmpdir, "blanks.gmt")
        out = os.path.join(self.tmpdir, "blanks.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATH1\tdesc1\tA\tB\n\n\nPATH2\tdesc2\tC\n")
        result = convert_gmt_to_engine_format(gmt, out)
        self.assertEqual(result["gene_sets_converted"], 2)

    def test_empty_description_handled(self):
        """GMT with empty description field (consecutive tabs) converts correctly."""
        gmt = os.path.join(self.tmpdir, "empty_desc.gmt")
        out = os.path.join(self.tmpdir, "empty_desc.engine.txt")
        with open(gmt, "w") as f:
            f.write("PATH1\t\tGENE1\tGENE2\n")
        result = convert_gmt_to_engine_format(gmt, out)
        self.assertEqual(result["gene_sets_converted"], 1)
        with open(out) as f:
            line = f.read().strip()
        self.assertEqual(line, "PATH1\tGENE1\tGENE2")


# ── Gene-Set File Validation ──

class TestGeneSetFileValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_valid_file(self):
        path = os.path.join(self.tmpdir, "valid.txt")
        with open(path, "w") as f:
            f.write("PATHWAY1\tGENE1\tGENE2\n")
        result = validate_gene_set_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_nonexistent_file(self):
        result = validate_gene_set_file("/does/not/exist")
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_empty_file(self):
        path = os.path.join(self.tmpdir, "empty.txt")
        with open(path, "w") as f:
            pass
        result = validate_gene_set_file(path)
        self.assertEqual(result.status, ValidationStatus.FAIL)


# ── Legal Combinations ──

class TestLegalCombinations(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a valid engine-format gene-set file
        self.valid_file = os.path.join(self.tmpdir, "custom.txt")
        with open(self.valid_file, "w") as f:
            f.write("MY_PATHWAY\tGENE1\tGENE2\tGENE3\n")
        self.outdir = os.path.join(self.tmpdir, "output")
        os.makedirs(self.outdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_default_no_custom_passes(self):
        """gene_sets=default, no custom files → PASS, empty paths."""
        r = prepare_custom_gene_sets([], None, None, "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)
        self.assertEqual(r.normalized_paths, [])

    def test_mouse_only_no_custom_passes(self):
        r = prepare_custom_gene_sets([], None, None, "mouse-only", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)

    def test_msigdb_only_no_custom_passes(self):
        r = prepare_custom_gene_sets([], None, None, "msigdb-only", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)

    def test_custom_profile_no_files_fails(self):
        """gene_sets=custom, no files → FAIL."""
        r = prepare_custom_gene_sets([], None, None, "custom", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)
        self.assertTrue(any("requires" in i for i in r.issues))

    def test_files_without_format_fails(self):
        """custom files provided, no format → FAIL."""
        r = prepare_custom_gene_sets([self.valid_file], None, "replace",
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)
        self.assertTrue(any("format" in i.lower() for i in r.issues))

    def test_files_without_action_fails(self):
        """custom files provided, no action → FAIL."""
        r = prepare_custom_gene_sets([self.valid_file], "engine", None,
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)
        self.assertTrue(any("action" in i.lower() for i in r.issues))

    def test_action_without_files_fails(self):
        """action=supplement, no files → FAIL."""
        r = prepare_custom_gene_sets([], None, "supplement",
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)

    def test_format_without_files_fails(self):
        """format=engine, no files → FAIL."""
        r = prepare_custom_gene_sets([], "engine", None,
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)

    def test_custom_replace_valid(self):
        """gene_sets=custom + files + format + action=replace → PASS."""
        r = prepare_custom_gene_sets([self.valid_file], "engine", "replace",
                                      "custom", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)
        self.assertEqual(len(r.normalized_paths), 1)
        self.assertEqual(r.gene_sets_copied, 1)

    def test_custom_supplement_valid(self):
        """gene_sets=custom + files + format + action=supplement → PASS."""
        r = prepare_custom_gene_sets([self.valid_file], "engine", "supplement",
                                      "custom", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)
        self.assertEqual(len(r.normalized_paths), 1)

    def test_default_with_supplement(self):
        """gene_sets=default + files + format + action=supplement → PASS."""
        r = prepare_custom_gene_sets([self.valid_file], "engine", "supplement",
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)

    def test_default_with_replace(self):
        """gene_sets=default + files + format + action=replace → PASS."""
        r = prepare_custom_gene_sets([self.valid_file], "engine", "replace",
                                      "default", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)

    def test_nonexistent_custom_file_fails(self):
        """Custom file that does not exist → FAIL."""
        r = prepare_custom_gene_sets(["/does/not/exist.txt"], "engine", "replace",
                                      "custom", self.outdir)
        self.assertEqual(r.status, ValidationStatus.FAIL)

    def test_gmt_conversion_in_prepare(self):
        """GMT file is converted to engine format during preparation."""
        gmt_file = os.path.join(self.tmpdir, "pathways.gmt")
        with open(gmt_file, "w") as f:
            f.write("PATH1\tdesc\tA\tB\n")
        r = prepare_custom_gene_sets([gmt_file], "gmt", "replace",
                                      "custom", self.outdir)
        self.assertEqual(r.status, ValidationStatus.PASS)
        self.assertEqual(r.gene_sets_converted, 1)
        # Verify the normalized file exists and has engine format
        self.assertEqual(len(r.normalized_paths), 1)
        with open(r.normalized_paths[0]) as f:
            line = f.read().strip()
        self.assertEqual(line, "PATH1\tA\tB")


# ── Resolve Gene Set Paths with Custom ──

class TestResolveGeneSetsWithCustom(unittest.TestCase):
    def setUp(self):
        self.base_dir = REPO_ROOT
        if not os.path.exists(os.path.join(REPO_ROOT, "data",
                                            "gene_set_list_mouse_2024.txt")):
            self.skipTest("Reference files not available")

    def test_default_without_custom_unchanged(self):
        """No custom paths → same as Phase 1."""
        paths = resolve_gene_set_paths("default", self.base_dir)
        self.assertEqual(len(paths), 2)
        self.assertIn("mouse_2024", paths[0])
        self.assertIn("msigdb_nohp", paths[1])

    def test_mouse_only_profile(self):
        """mouse-only → 1 path."""
        paths = resolve_gene_set_paths("mouse-only", self.base_dir)
        self.assertEqual(len(paths), 1)
        self.assertIn("mouse_2024", paths[0])

    def test_msigdb_only_profile(self):
        """msigdb-only → 1 path."""
        paths = resolve_gene_set_paths("msigdb-only", self.base_dir)
        self.assertEqual(len(paths), 1)
        self.assertIn("msigdb_nohp", paths[0])

    def test_custom_replace_returns_only_custom(self):
        """action=replace → only custom paths."""
        custom = ["/tmp/custom1.txt", "/tmp/custom2.txt"]
        paths = resolve_gene_set_paths("default", self.base_dir,
                                        custom_normalized_paths=custom,
                                        custom_action="replace")
        self.assertEqual(paths, custom)

    def test_custom_supplement_appends(self):
        """action=supplement → built-in + custom."""
        custom = ["/tmp/custom1.txt"]
        paths = resolve_gene_set_paths("default", self.base_dir,
                                        custom_normalized_paths=custom,
                                        custom_action="supplement")
        self.assertEqual(len(paths), 3)  # 2 default + 1 custom
        self.assertEqual(paths[-1], "/tmp/custom1.txt")

    def test_mouse_only_supplement_appends(self):
        """mouse-only + supplement → mouse + custom."""
        custom = ["/tmp/custom1.txt"]
        paths = resolve_gene_set_paths("mouse-only", self.base_dir,
                                        custom_normalized_paths=custom,
                                        custom_action="supplement")
        self.assertEqual(len(paths), 2)  # 1 mouse + 1 custom
        self.assertIn("mouse_2024", paths[0])
        self.assertEqual(paths[1], "/tmp/custom1.txt")

    def test_custom_profile_replace(self):
        """gene_sets=custom + replace → only custom."""
        custom = ["/tmp/custom1.txt"]
        paths = resolve_gene_set_paths("custom", self.base_dir,
                                        custom_normalized_paths=custom,
                                        custom_action="replace")
        self.assertEqual(paths, custom)

    def test_custom_profile_supplement(self):
        """gene_sets=custom + supplement → default base + custom."""
        custom = ["/tmp/custom1.txt"]
        paths = resolve_gene_set_paths("custom", self.base_dir,
                                        custom_normalized_paths=custom,
                                        custom_action="supplement")
        self.assertEqual(len(paths), 3)  # 2 default + 1 custom


# ── Background Validation ──

class TestBackgroundValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        if not os.path.exists(GENE_MAP):
            self.skipTest("Reference files not available")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_all_input_in_background(self):
        """All input genes present in background → PASS."""
        from pigean.validation import validate_background
        bg_path = os.path.join(self.tmpdir, "bg.txt")
        input_path = os.path.join(self.tmpdir, "input.txt")
        with open(bg_path, "w") as f:
            f.write("LEP\nADIPOQ\nBRD2\n")
        with open(input_path, "w") as f:
            f.write("LEP\nADIPOQ\n")
        result = validate_background(bg_path, input_path, GENE_MAP)
        self.assertEqual(result.status, ValidationStatus.PASS)
        self.assertEqual(result.input_genes_missing_from_background, [])

    def test_some_input_missing_from_background(self):
        """Some input genes not in background → PASS_WITH_WARNINGS."""
        from pigean.validation import validate_background
        bg_path = os.path.join(self.tmpdir, "bg.txt")
        input_path = os.path.join(self.tmpdir, "input.txt")
        with open(bg_path, "w") as f:
            f.write("LEP\nBRD2\n")
        with open(input_path, "w") as f:
            f.write("LEP\nADIPOQ\n")
        result = validate_background(bg_path, input_path, GENE_MAP)
        self.assertEqual(result.status, ValidationStatus.PASS_WITH_WARNINGS)
        self.assertIn("ADIPOQ", result.input_genes_missing_from_background)

    def test_no_recognized_background_genes(self):
        """Zero recognized background genes → FAIL."""
        from pigean.validation import validate_background
        bg_path = os.path.join(self.tmpdir, "bg.txt")
        input_path = os.path.join(self.tmpdir, "input.txt")
        with open(bg_path, "w") as f:
            f.write("FAKEGENE1\nFAKEGENE2\n")
        with open(input_path, "w") as f:
            f.write("LEP\n")
        result = validate_background(bg_path, input_path, GENE_MAP)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_background_larger_than_input(self):
        """Background has many more genes than input → PASS."""
        from pigean.validation import validate_background
        bg_path = os.path.join(self.tmpdir, "bg.txt")
        input_path = os.path.join(self.tmpdir, "input.txt")
        with open(bg_path, "w") as f:
            f.write("LEP\nADIPOQ\nBRD2\nSOX2\nLITAF\nSLCO1B1\n")
        with open(input_path, "w") as f:
            f.write("LEP\n")
        result = validate_background(bg_path, input_path, GENE_MAP)
        self.assertEqual(result.status, ValidationStatus.PASS)
        self.assertEqual(result.total_count, 6)
        self.assertGreater(result.mapped_count, 0)


if __name__ == "__main__":
    unittest.main()
