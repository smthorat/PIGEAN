"""Tests for Phase 3 evidence adapters."""
import sys, os, unittest, tempfile, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.adapters import (
    ValidationStatus, get_adapter, SUPPORTED_ANALYSIS_TYPES
)
from pigean.adapters.positive_controls import PositiveControlsAdapter
from pigean.adapters.bayes_factor import BayesFactorAdapter
from pigean.adapters.zscore import ZScoreAdapter
from pigean.adapters.percentile import PercentileAdapter
from pigean.adapters.exome import ExomeAdapter

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "evidence")


class TestAdapterRegistry(unittest.TestCase):
    def test_positive_controls_returns_correct_adapter(self):
        adapter = get_adapter("positive-controls")
        self.assertIsInstance(adapter, PositiveControlsAdapter)

    def test_bayes_factor_returns_correct_adapter(self):
        adapter = get_adapter("gene-bayes-factor")
        self.assertIsInstance(adapter, BayesFactorAdapter)

    def test_zscore_returns_correct_adapter(self):
        adapter = get_adapter("gene-z-score")
        self.assertIsInstance(adapter, ZScoreAdapter)

    def test_percentile_returns_correct_adapter(self):
        adapter = get_adapter("gene-percentile")
        self.assertIsInstance(adapter, PercentileAdapter)

    def test_exome_returns_correct_adapter(self):
        adapter = get_adapter("exome")
        self.assertIsInstance(adapter, ExomeAdapter)

    def test_unknown_analysis_raises(self):
        with self.assertRaises(ValueError):
            get_adapter("unknown-type")

    def test_all_supported_types_have_adapters(self):
        for t in SUPPORTED_ANALYSIS_TYPES:
            adapter = get_adapter(t)
            self.assertIsNotNone(adapter)


class TestBayesFactorAdapter(unittest.TestCase):
    def setUp(self):
        self.adapter = BayesFactorAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_validate_valid_file(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_valid.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_validate_missing_column(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_invalid_nocol.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("log_bf" in i for i in result.issues))

    def test_validate_custom_columns(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_invalid_nocol.tsv")
        params = {"gene_column": "Gene", "score_column": "score"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_validate_nonexistent(self):
        result = self.adapter.validate_file("/does/not/exist.tsv")
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_validate_empty(self):
        path = os.path.join(self.tmpdir, "empty.tsv")
        open(path, "w").close()
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_normalize_preserves_original(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_valid.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertTrue(os.path.exists(result.original_path))
        self.assertTrue(os.path.exists(result.normalized_path))

    def test_normalize_counts(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_valid.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertEqual(result.original_rows, 6)
        self.assertEqual(result.normalized_rows, 6)
        self.assertEqual(result.removed_rows, 0)

    def test_normalize_removes_na(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_with_na.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertEqual(result.original_rows, 4)
        self.assertEqual(result.normalized_rows, 3)
        self.assertEqual(result.removed_rows, 1)

    def test_normalize_has_transformations(self):
        path = os.path.join(FIXTURES_DIR, "bayes_factor_valid.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertTrue(len(result.transformations) > 0)
        self.assertTrue(any("natural-log" in t for t in result.transformations))

    def test_build_engine_args(self):
        args = self.adapter.build_engine_args("/tmp/bf.tsv")
        self.assertEqual(args, ["--gene-bfs-in", "/tmp/bf.tsv"])

    def test_engine_flag_attribute(self):
        self.assertEqual(self.adapter.ENGINE_FLAG, "--gene-bfs-in")


class TestZScoreAdapter(unittest.TestCase):
    def setUp(self):
        self.adapter = ZScoreAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_validate_requires_gene_column(self):
        path = os.path.join(FIXTURES_DIR, "zscore_valid.tsv")
        result = self.adapter.validate_file(path, {"score_column": "zscore"})
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("gene-column" in i.lower() for i in result.issues))

    def test_validate_requires_score_column(self):
        path = os.path.join(FIXTURES_DIR, "zscore_valid.tsv")
        result = self.adapter.validate_file(path, {"gene_column": "Gene"})
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("score-column" in i.lower() for i in result.issues))

    def test_validate_valid(self):
        path = os.path.join(FIXTURES_DIR, "zscore_valid.tsv")
        params = {"gene_column": "Gene", "score_column": "zscore"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_validate_wrong_column_name(self):
        path = os.path.join(FIXTURES_DIR, "zscore_valid.tsv")
        params = {"gene_column": "Gene", "score_column": "NONEXISTENT"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_normalize(self):
        path = os.path.join(FIXTURES_DIR, "zscore_valid.tsv")
        params = {"gene_column": "Gene", "score_column": "zscore"}
        result = self.adapter.normalize(path, self.tmpdir, params)
        self.assertEqual(result.normalized_rows, 4)
        self.assertTrue(os.path.exists(result.normalized_path))

    def test_build_engine_args(self):
        params = {"gene_column": "Gene", "score_column": "zscore"}
        args = self.adapter.build_engine_args("/tmp/zs.tsv", params)
        self.assertIn("--gene-zs-in", args)
        self.assertIn("--gene-zs-id-col", args)
        idx = args.index("--gene-zs-id-col")
        self.assertEqual(args[idx + 1], "Gene")


class TestPercentileAdapter(unittest.TestCase):
    def setUp(self):
        self.adapter = PercentileAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_validate_requires_gene_column(self):
        path = os.path.join(FIXTURES_DIR, "percentile_valid.tsv")
        result = self.adapter.validate_file(path, {"score_column": "rank_score"})
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_validate_valid(self):
        path = os.path.join(FIXTURES_DIR, "percentile_valid.tsv")
        params = {"gene_column": "Gene", "score_column": "rank_score"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_normalize(self):
        path = os.path.join(FIXTURES_DIR, "percentile_valid.tsv")
        params = {"gene_column": "Gene", "score_column": "rank_score"}
        result = self.adapter.normalize(path, self.tmpdir, params)
        self.assertEqual(result.normalized_rows, 5)

    def test_build_engine_args_default(self):
        params = {"gene_column": "Gene", "score_column": "pct"}
        args = self.adapter.build_engine_args("/tmp/p.tsv", params)
        self.assertNotIn("--gene-percentiles-higher-is-better", args)

    def test_build_engine_args_higher_is_better(self):
        params = {"gene_column": "Gene", "score_column": "pct", "higher_is_better": True}
        args = self.adapter.build_engine_args("/tmp/p.tsv", params)
        self.assertIn("--gene-percentiles-higher-is-better", args)


class TestExomeAdapter(unittest.TestCase):
    def setUp(self):
        self.adapter = ExomeAdapter()
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_validate_valid(self):
        path = os.path.join(FIXTURES_DIR, "exome_valid.tsv")
        result = self.adapter.validate_file(path)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_validate_nonexistent(self):
        result = self.adapter.validate_file("/does/not/exist.tsv")
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_validate_wrong_specified_column(self):
        path = os.path.join(FIXTURES_DIR, "exome_valid.tsv")
        params = {"exomes_gene_col": "NONEXISTENT"}
        result = self.adapter.validate_file(path, params)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_normalize(self):
        path = os.path.join(FIXTURES_DIR, "exome_valid.tsv")
        result = self.adapter.normalize(path, self.tmpdir)
        self.assertEqual(result.normalized_rows, 5)
        self.assertTrue(os.path.exists(result.original_path))
        self.assertTrue(os.path.exists(result.normalized_path))

    def test_build_engine_args_minimal(self):
        args = self.adapter.build_engine_args("/tmp/ex.tsv")
        self.assertEqual(args, ["--exomes-in", "/tmp/ex.tsv"])

    def test_build_engine_args_with_overrides(self):
        params = {
            "exomes_gene_col": "GENE",
            "exomes_p_col": "PVAL",
            "exomes_n": 50000,
        }
        args = self.adapter.build_engine_args("/tmp/ex.tsv", params)
        self.assertIn("--exomes-gene-col", args)
        self.assertIn("--exomes-n", args)
        idx = args.index("--exomes-n")
        self.assertEqual(args[idx + 1], "50000")


class TestPositiveControlsBackwardCompat(unittest.TestCase):
    """Verify positive-controls adapter still works identically."""

    def setUp(self):
        self.adapter = PositiveControlsAdapter()

    def test_is_base_evidence_adapter(self):
        from pigean.adapters.base import BaseEvidenceAdapter
        self.assertIsInstance(self.adapter, BaseEvidenceAdapter)

    def test_analysis_type(self):
        self.assertEqual(self.adapter.ANALYSIS_TYPE, "positive-controls")

    def test_engine_flag(self):
        self.assertEqual(self.adapter.ENGINE_FLAG, "--positive-controls-in")

    def test_build_engine_args(self):
        args = self.adapter.build_engine_args("/tmp/genes.txt")
        self.assertEqual(args, ["--positive-controls-in", "/tmp/genes.txt"])


if __name__ == "__main__":
    unittest.main()
