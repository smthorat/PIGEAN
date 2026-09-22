"""Tests for pigean/engine.py — verify the constructed command matches Phase 0."""
import sys, os, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.config import DEFAULTS
from pigean.engine import build_priors_command

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class TestEngineCommand(unittest.TestCase):
    def setUp(self):
        self.config = dict(DEFAULTS)
        self.base_dir = REPO_ROOT
        self.output_dir = "/tmp/test_output"
        self.normalized_path = "/tmp/test_output/normalized/positive_controls.txt"
        self.evidence_args = ["--positive-controls-in", self.normalized_path]
        self.reference_paths = {
            "gene_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.loc"),
            "tss_loc": os.path.join(REPO_ROOT, "data", "refGene_hg19_TSS.subset.loc"),
            "exons_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.exons.loc"),
        }
        self.gene_set_paths = [
            os.path.join(REPO_ROOT, "data", "gene_set_list_mouse_2024.txt"),
            os.path.join(REPO_ROOT, "data", "gene_set_list_msigdb_nohp.txt"),
        ]
        self.gene_map_path = os.path.join(REPO_ROOT, "data", "portal_gencode.gene.map")

    def _build(self):
        return build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir, self.base_dir
        )

    def test_starts_with_python_priors_gibbs(self):
        cmd = self._build()
        self.assertEqual(cmd[0], "python3")
        self.assertEqual(cmd[1], "-u")
        self.assertTrue(cmd[2].endswith("priors.py"))
        self.assertEqual(cmd[3], "gibbs")

    def test_gene_set_order_matches_phase0(self):
        cmd = self._build()
        x_in_indices = [i for i, v in enumerate(cmd) if v == "--X-in"]
        self.assertEqual(len(x_in_indices), 2)
        first_file = cmd[x_in_indices[0] + 1]
        second_file = cmd[x_in_indices[1] + 1]
        self.assertIn("mouse_2024", first_file)
        self.assertIn("msigdb_nohp", second_file)

    def test_uses_full_flag_names(self):
        cmd = self._build()
        self.assertIn("--gene-gene-set-stats-out", cmd)

    def test_correct_gene_map(self):
        cmd = self._build()
        idx = cmd.index("--gene-map-in")
        self.assertIn("portal_gencode.gene.map", cmd[idx + 1])

    def test_correct_references(self):
        cmd = self._build()
        idx = cmd.index("--gene-loc-file")
        self.assertIn("NCBI37.3.plink.gene.loc", cmd[idx + 1])
        idx = cmd.index("--gene-loc-file-huge")
        self.assertIn("refGene_hg19_TSS.subset.loc", cmd[idx + 1])
        idx = cmd.index("--exons-loc-file-huge")
        self.assertIn("NCBI37.3.plink.gene.exons.loc", cmd[idx + 1])

    def test_filter_values_match_phase0(self):
        cmd = self._build()
        idx = cmd.index("--gene-filter-value")
        self.assertEqual(cmd[idx + 1], "1")
        idx = cmd.index("--gene-set-filter-value")
        self.assertEqual(cmd[idx + 1], "0.01")

    def test_max_gene_sets(self):
        cmd = self._build()
        idx = cmd.index("--max-num-gene-sets")
        self.assertEqual(cmd[idx + 1], "5000")

    def test_debug_level(self):
        cmd = self._build()
        idx = cmd.index("--debug-level")
        self.assertEqual(cmd[idx + 1], "3")

    def test_positive_controls_input(self):
        cmd = self._build()
        idx = cmd.index("--positive-controls-in")
        self.assertEqual(cmd[idx + 1], self.normalized_path)

    def test_output_files_standardized(self):
        cmd = self._build()
        idx = cmd.index("--gene-stats-out")
        self.assertTrue(cmd[idx + 1].endswith("gs.out"))
        idx = cmd.index("--gene-set-stats-out")
        self.assertTrue(cmd[idx + 1].endswith("gss.out"))
        idx = cmd.index("--gene-gene-set-stats-out")
        self.assertTrue(cmd[idx + 1].endswith("ggss.out"))
        idx = cmd.index("--params-out")
        self.assertTrue(cmd[idx + 1].endswith("p.out"))

    def test_no_extra_scientific_flags(self):
        cmd = self._build()
        self.assertNotIn("--max-num-iter", cmd)
        self.assertNotIn("--num-chains", cmd)
        self.assertNotIn("--sparse-solution", cmd)

    def test_no_background_flag_by_default(self):
        cmd = self._build()
        self.assertNotIn("--positive-controls-all-in", cmd)


class TestBackgroundCommand(unittest.TestCase):
    def setUp(self):
        self.config = dict(DEFAULTS)
        self.base_dir = REPO_ROOT
        self.output_dir = "/tmp/test_output"
        self.normalized_path = "/tmp/test_output/normalized/positive_controls.txt"
        self.evidence_args = ["--positive-controls-in", self.normalized_path]
        self.reference_paths = {
            "gene_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.loc"),
            "tss_loc": os.path.join(REPO_ROOT, "data", "refGene_hg19_TSS.subset.loc"),
            "exons_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.exons.loc"),
        }
        self.gene_set_paths = [
            os.path.join(REPO_ROOT, "data", "gene_set_list_mouse_2024.txt"),
            os.path.join(REPO_ROOT, "data", "gene_set_list_msigdb_nohp.txt"),
        ]
        self.gene_map_path = os.path.join(REPO_ROOT, "data", "portal_gencode.gene.map")

    def test_with_background_adds_flag(self):
        bg_path = "/tmp/normalized/background.txt"
        cmd = build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir, background_normalized_path=bg_path
        )
        self.assertIn("--positive-controls-all-in", cmd)
        idx = cmd.index("--positive-controls-all-in")
        self.assertEqual(cmd[idx + 1], bg_path)

    def test_background_does_not_affect_other_flags(self):
        cmd_no_bg = build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir
        )
        cmd_with_bg = build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir, background_normalized_path="/tmp/bg.txt"
        )
        bg_idx = cmd_with_bg.index("--positive-controls-all-in")
        cmd_without_bg_tokens = cmd_with_bg[:bg_idx] + cmd_with_bg[bg_idx + 2:]
        self.assertEqual(cmd_no_bg, cmd_without_bg_tokens)

    def test_explicit_none_background_same_as_omitted(self):
        cmd_default = build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir
        )
        cmd_none = build_priors_command(
            self.config, self.evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir, background_normalized_path=None
        )
        self.assertEqual(cmd_default, cmd_none)


class TestEvidenceAdapterEngineArgs(unittest.TestCase):
    """Test that each adapter produces the correct engine arguments."""

    def setUp(self):
        self.config = dict(DEFAULTS)
        self.base_dir = REPO_ROOT
        self.output_dir = "/tmp/test_output"
        self.reference_paths = {
            "gene_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.loc"),
            "tss_loc": os.path.join(REPO_ROOT, "data", "refGene_hg19_TSS.subset.loc"),
            "exons_loc": os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.exons.loc"),
        }
        self.gene_set_paths = [
            os.path.join(REPO_ROOT, "data", "gene_set_list_mouse_2024.txt"),
        ]
        self.gene_map_path = os.path.join(REPO_ROOT, "data", "portal_gencode.gene.map")

    def test_bayes_factor_args(self):
        from pigean.adapters.bayes_factor import BayesFactorAdapter
        adapter = BayesFactorAdapter()
        args = adapter.build_engine_args("/tmp/bf.tsv")
        self.assertEqual(args, ["--gene-bfs-in", "/tmp/bf.tsv"])

    def test_zscore_args(self):
        from pigean.adapters.zscore import ZScoreAdapter
        adapter = ZScoreAdapter()
        params = {"gene_column": "Gene", "score_column": "z"}
        args = adapter.build_engine_args("/tmp/zs.tsv", params)
        self.assertIn("--gene-zs-in", args)
        self.assertIn("--gene-zs-id-col", args)
        self.assertIn("--gene-zs-value-col", args)
        idx = args.index("--gene-zs-id-col")
        self.assertEqual(args[idx + 1], "Gene")
        idx = args.index("--gene-zs-value-col")
        self.assertEqual(args[idx + 1], "z")

    def test_percentile_args(self):
        from pigean.adapters.percentile import PercentileAdapter
        adapter = PercentileAdapter()
        params = {"gene_column": "Gene", "score_column": "pct"}
        args = adapter.build_engine_args("/tmp/pct.tsv", params)
        self.assertIn("--gene-percentiles-in", args)
        self.assertIn("--gene-percentiles-id-col", args)
        self.assertIn("--gene-percentiles-value-col", args)

    def test_percentile_higher_is_better(self):
        from pigean.adapters.percentile import PercentileAdapter
        adapter = PercentileAdapter()
        params = {"gene_column": "G", "score_column": "V", "higher_is_better": True}
        args = adapter.build_engine_args("/tmp/pct.tsv", params)
        self.assertIn("--gene-percentiles-higher-is-better", args)

    def test_percentile_default_lower_is_better(self):
        from pigean.adapters.percentile import PercentileAdapter
        adapter = PercentileAdapter()
        params = {"gene_column": "G", "score_column": "V"}
        args = adapter.build_engine_args("/tmp/pct.tsv", params)
        self.assertNotIn("--gene-percentiles-higher-is-better", args)

    def test_exome_args_minimal(self):
        from pigean.adapters.exome import ExomeAdapter
        adapter = ExomeAdapter()
        args = adapter.build_engine_args("/tmp/ex.tsv")
        self.assertEqual(args, ["--exomes-in", "/tmp/ex.tsv"])

    def test_exome_args_with_columns(self):
        from pigean.adapters.exome import ExomeAdapter
        adapter = ExomeAdapter()
        params = {"exomes_gene_col": "GENE", "exomes_p_col": "PVALUE"}
        args = adapter.build_engine_args("/tmp/ex.tsv", params)
        self.assertIn("--exomes-gene-col", args)
        self.assertIn("--exomes-p-col", args)

    def test_bf_args_in_full_command(self):
        evidence_args = ["--gene-bfs-in", "/tmp/bf.tsv"]
        cmd = build_priors_command(
            self.config, evidence_args, self.reference_paths,
            self.gene_set_paths, self.gene_map_path, self.output_dir,
            self.base_dir
        )
        self.assertIn("--gene-bfs-in", cmd)
        self.assertNotIn("--positive-controls-in", cmd)


if __name__ == "__main__":
    unittest.main()
