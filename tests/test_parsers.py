"""Tests for pigean/parsers.py — Phase 5 output parsers."""

import os
import tempfile
import unittest

from pigean.parsers import (
    parse_gene_stats,
    parse_gene_set_stats,
    parse_gene_gene_set_stats,
    parse_params,
)

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "outputs")


class TestParseGeneStatsPositiveControls(unittest.TestCase):
    """Tests for parse_gene_stats with positive-controls output."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "gs_positive_controls.out")
        self.result = parse_gene_stats(self.path)

    def test_returns_dict(self):
        self.assertIsNotNone(self.result)
        self.assertIsInstance(self.result, dict)

    def test_row_count(self):
        self.assertEqual(self.result["row_count"], 5)

    def test_columns_detected(self):
        expected = ["Gene", "prior", "prior_adj", "combined", "combined_adj",
                    "combined_D", "positive_control", "log_bf", "N",
                    "Chrom", "Start", "End"]
        self.assertEqual(self.result["columns"], expected)

    def test_has_positive_control(self):
        self.assertTrue(self.result["has_positive_control"])

    def test_sorted_by_combined_D_descending(self):
        genes = self.result["genes"]
        vals = [g["combined_D"] for g in genes]
        self.assertEqual(vals, sorted(vals, reverse=True))

    def test_first_gene(self):
        top = self.result["genes"][0]
        self.assertEqual(top["Gene"], "GENE_A")
        self.assertAlmostEqual(top["combined_D"], 0.987)
        self.assertAlmostEqual(top["prior"], 2.09)
        self.assertAlmostEqual(top["log_bf"], 5.89)
        self.assertEqual(top["N"], 930)

    def test_positive_control_values(self):
        genes = self.result["genes"]
        pc_vals = {g["Gene"]: g["positive_control"] for g in genes}
        self.assertAlmostEqual(pc_vals["GENE_A"], 5.89)
        self.assertAlmostEqual(pc_vals["GENE_D"], 0.0)

    def test_float_conversion(self):
        gene = self.result["genes"][0]
        self.assertIsInstance(gene["prior"], float)
        self.assertIsInstance(gene["combined_D"], float)
        self.assertIsInstance(gene["log_bf"], float)

    def test_int_conversion(self):
        gene = self.result["genes"][0]
        self.assertIsInstance(gene["N"], int)
        self.assertIsInstance(gene["Start"], int)
        self.assertIsInstance(gene["End"], int)

    def test_string_columns(self):
        gene = self.result["genes"][0]
        self.assertIsInstance(gene["Gene"], str)
        self.assertIsInstance(gene["Chrom"], str)

    def test_location_data(self):
        gene = self.result["genes"][0]
        self.assertEqual(gene["Chrom"], "7")
        self.assertEqual(gene["Start"], 127881241)
        self.assertEqual(gene["End"], 127897682)


class TestParseGeneStatsBFMode(unittest.TestCase):
    """Tests for parse_gene_stats with BF mode output (no positive_control)."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "gs_bf_mode.out")
        self.result = parse_gene_stats(self.path)

    def test_no_positive_control(self):
        self.assertFalse(self.result["has_positive_control"])

    def test_column_count(self):
        self.assertEqual(len(self.result["columns"]), 11)

    def test_row_count(self):
        self.assertEqual(self.result["row_count"], 5)


class TestParseGeneStatsMissingFile(unittest.TestCase):
    """Tests for parse_gene_stats with missing file."""

    def test_missing_file_returns_none(self):
        result = parse_gene_stats("/tmp/nonexistent_gs.out")
        self.assertIsNone(result)

    def test_empty_file_returns_none(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".out",
                                          delete=False) as f:
            f.write("")
            f.flush()
            try:
                result = parse_gene_stats(f.name)
                self.assertIsNone(result)
            finally:
                os.unlink(f.name)

    def test_header_only_returns_zero_rows(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".out",
                                          delete=False) as f:
            f.write("Gene\tprior\tcombined_D\n")
            f.flush()
            try:
                result = parse_gene_stats(f.name)
                self.assertIsNotNone(result)
                self.assertEqual(result["row_count"], 0)
            finally:
                os.unlink(f.name)


class TestParseGeneSetStats(unittest.TestCase):
    """Tests for parse_gene_set_stats."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "gss_sample.out")
        self.result = parse_gene_set_stats(self.path)

    def test_returns_dict(self):
        self.assertIsNotNone(self.result)

    def test_row_count(self):
        self.assertEqual(self.result["row_count"], 5)

    def test_sorted_by_abs_beta(self):
        gene_sets = self.result["gene_sets"]
        betas = [abs(gs["beta"]) for gs in gene_sets]
        self.assertEqual(betas, sorted(betas, reverse=True))

    def test_first_gene_set_highest_abs_beta(self):
        top = self.result["gene_sets"][0]
        # gs_pathway_alpha has beta=0.0219, highest absolute value
        self.assertEqual(top["Gene_Set"], "gs_pathway_alpha")
        self.assertAlmostEqual(top["beta"], 0.0219)

    def test_column_count(self):
        self.assertEqual(len(self.result["columns"]), 26)

    def test_verified_fields(self):
        gs = self.result["gene_sets"][0]
        self.assertIn("Gene_Set", gs)
        self.assertIn("label", gs)
        self.assertIn("N", gs)
        self.assertIn("beta", gs)
        self.assertIn("avg_postp", gs)
        self.assertIn("P", gs)

    def test_gene_set_name_is_string(self):
        gs = self.result["gene_sets"][0]
        self.assertIsInstance(gs["Gene_Set"], str)
        self.assertIsInstance(gs["label"], str)

    def test_missing_file(self):
        result = parse_gene_set_stats("/tmp/nonexistent_gss.out")
        self.assertIsNone(result)


class TestParseGeneGeneSetStats(unittest.TestCase):
    """Tests for parse_gene_gene_set_stats with bidirectional index."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "ggss_sample.out")
        self.result = parse_gene_gene_set_stats(self.path)

    def test_returns_dict(self):
        self.assertIsNotNone(self.result)

    def test_row_count(self):
        self.assertEqual(self.result["row_count"], 10)

    def test_by_gene_index(self):
        by_gene = self.result["by_gene"]
        self.assertIn("GENE_A", by_gene)
        # GENE_A appears in 3 gene sets
        self.assertEqual(len(by_gene["GENE_A"]), 3)

    def test_by_gene_set_index(self):
        by_gs = self.result["by_gene_set"]
        self.assertIn("gs_pathway_alpha", by_gs)
        # gs_pathway_alpha has 2 genes
        self.assertEqual(len(by_gs["gs_pathway_alpha"]), 2)

    def test_entry_fields(self):
        entry = self.result["entries"][0]
        self.assertIn("Gene", entry)
        self.assertIn("gene_set", entry)
        self.assertIn("beta", entry)
        self.assertIn("weight", entry)

    def test_all_genes_indexed(self):
        by_gene = self.result["by_gene"]
        self.assertEqual(set(by_gene.keys()),
                         {"GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"})

    def test_missing_file(self):
        result = parse_gene_gene_set_stats("/tmp/nonexistent_ggss.out")
        self.assertIsNone(result)


class TestParseParamsConverged(unittest.TestCase):
    """Tests for parse_params with converged output."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "p_converged.out")
        self.result = parse_params(self.path)

    def test_returns_dict(self):
        self.assertIsNotNone(self.result)

    def test_has_params_and_raw(self):
        self.assertIn("params", self.result)
        self.assertIn("raw_rows", self.result)

    def test_num_gibbs_iter(self):
        self.assertEqual(self.result["params"]["num_gibbs_iter"], 200)

    def test_num_restarts(self):
        self.assertEqual(self.result["params"]["num_gibbs_restarts"], 0)

    def test_num_chains(self):
        self.assertEqual(self.result["params"]["num_chains"], 10)

    def test_float_param(self):
        p = self.result["params"]["p"]
        self.assertIsInstance(p, float)
        self.assertAlmostEqual(p, 0.00342)

    def test_boolean_param(self):
        self.assertTrue(self.result["params"]["use_mean_betas"])
        self.assertTrue(self.result["params"]["sparse_solution"])

    def test_single_version_is_scalar(self):
        # Single-version params stored as scalars
        val = self.result["params"]["num_chains"]
        self.assertNotIsInstance(val, list)


class TestParseParamsMultiVersion(unittest.TestCase):
    """Tests for parse_params with multi-version parameters."""

    def setUp(self):
        self.path = os.path.join(FIXTURES, "p_hit_cap.out")
        self.result = parse_params(self.path)

    def test_multi_version_p_stored_as_list(self):
        p = self.result["params"]["p"]
        self.assertIsInstance(p, list)
        self.assertEqual(len(p), 2)
        self.assertAlmostEqual(p[0], 0.00342)
        self.assertAlmostEqual(p[1], 0.00189)

    def test_multi_version_sigma2_stored_as_list(self):
        s2 = self.result["params"]["sigma2"]
        self.assertIsInstance(s2, list)
        self.assertEqual(len(s2), 2)

    def test_hit_cap_iter(self):
        self.assertEqual(self.result["params"]["num_gibbs_iter"], 499)

    def test_missing_file(self):
        result = parse_params("/tmp/nonexistent_p.out")
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
