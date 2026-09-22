"""Tests for pigean/interpretation.py — Phase 5 interpretation module."""

import os
import unittest

from pigean.parsers import (
    parse_gene_stats,
    parse_gene_set_stats,
    parse_gene_gene_set_stats,
    parse_params,
)
from pigean.convergence import assess_convergence
from pigean.interpretation import (
    FORBIDDEN_TERMS,
    INTERPRETATION_CONTRACTS,
    interpret_results,
)

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "outputs")


def _get_full_interpretation():
    """Load fixtures and produce a full interpretation result."""
    gs = parse_gene_stats(os.path.join(FIXTURES, "gs_positive_controls.out"))
    gss = parse_gene_set_stats(os.path.join(FIXTURES, "gss_sample.out"))
    ggss = parse_gene_gene_set_stats(os.path.join(FIXTURES, "ggss_sample.out"))
    params = parse_params(os.path.join(FIXTURES, "p_converged.out"))
    convergence = assess_convergence(params)
    return interpret_results(gs, gss, ggss, params, convergence,
                             "positive-controls")


class TestInterpretResultsStructure(unittest.TestCase):
    """Test the overall structure of interpret_results output."""

    def setUp(self):
        self.result = _get_full_interpretation()

    def test_has_all_keys(self):
        expected = {"gene_summary", "gene_set_summary", "gene_pathway_links",
                    "convergence_note", "interpretation_contracts",
                    "caveats", "analysis_type"}
        self.assertTrue(expected.issubset(set(self.result.keys())))

    def test_analysis_type(self):
        self.assertEqual(self.result["analysis_type"], "positive-controls")

    def test_gene_summary_not_none(self):
        self.assertIsNotNone(self.result["gene_summary"])

    def test_gene_set_summary_not_none(self):
        self.assertIsNotNone(self.result["gene_set_summary"])

    def test_gene_pathway_links_not_none(self):
        self.assertIsNotNone(self.result["gene_pathway_links"])

    def test_caveats_not_empty(self):
        self.assertTrue(len(self.result["caveats"]) > 0)


class TestInputCandidateSeparation(unittest.TestCase):
    """Test separation of input genes from candidate genes."""

    def setUp(self):
        self.result = _get_full_interpretation()
        self.gs = self.result["gene_summary"]

    def test_input_genes_identified(self):
        input_genes = self.gs["input_genes"]
        # Fixture has 3 input genes: GENE_A, GENE_B, GENE_C (pc=5.89)
        self.assertEqual(len(input_genes), 3)

    def test_candidate_genes_identified(self):
        candidate_genes = self.gs["candidate_genes"]
        # Fixture has 2 candidates: GENE_D, GENE_E (pc=0)
        self.assertEqual(len(candidate_genes), 2)

    def test_input_gene_status(self):
        for g in self.gs["input_genes"]:
            self.assertEqual(g["status"], "INPUT")

    def test_candidate_gene_status(self):
        for g in self.gs["candidate_genes"]:
            self.assertEqual(g["status"], "CANDIDATE")

    def test_input_gene_names(self):
        names = {g["gene"] for g in self.gs["input_genes"]}
        self.assertEqual(names, {"GENE_A", "GENE_B", "GENE_C"})

    def test_candidate_gene_names(self):
        names = {g["gene"] for g in self.gs["candidate_genes"]}
        self.assertEqual(names, {"GENE_D", "GENE_E"})


class TestBFModeNoCandidateSplit(unittest.TestCase):
    """BF mode has no positive_control column — all genes are candidates."""

    def setUp(self):
        gs = parse_gene_stats(os.path.join(FIXTURES, "gs_bf_mode.out"))
        gss = parse_gene_set_stats(os.path.join(FIXTURES, "gss_sample.out"))
        ggss = parse_gene_gene_set_stats(
            os.path.join(FIXTURES, "ggss_sample.out"))
        params = parse_params(os.path.join(FIXTURES, "p_converged.out"))
        convergence = assess_convergence(params)
        self.result = interpret_results(gs, gss, ggss, params, convergence,
                                        "gene-bayes-factor")

    def test_no_input_genes(self):
        self.assertEqual(len(self.result["gene_summary"]["input_genes"]), 0)

    def test_all_candidates(self):
        self.assertEqual(
            len(self.result["gene_summary"]["candidate_genes"]), 5)


class TestGeneSummaryRanking(unittest.TestCase):
    """Test that genes are ranked by combined_D descending."""

    def setUp(self):
        self.result = _get_full_interpretation()

    def test_input_genes_ranked(self):
        input_genes = self.result["gene_summary"]["input_genes"]
        vals = [g["combined_D"] for g in input_genes]
        self.assertEqual(vals, sorted(vals, reverse=True))

    def test_candidate_genes_ranked(self):
        candidate_genes = self.result["gene_summary"]["candidate_genes"]
        vals = [g["combined_D"] for g in candidate_genes]
        self.assertEqual(vals, sorted(vals, reverse=True))

    def test_rank_numbers_sequential(self):
        input_genes = self.result["gene_summary"]["input_genes"]
        ranks = [g["rank"] for g in input_genes]
        self.assertEqual(ranks, list(range(1, len(ranks) + 1)))


class TestGeneEntryFields(unittest.TestCase):
    """Test that gene entries have the correct VERIFIED fields."""

    def setUp(self):
        self.result = _get_full_interpretation()
        self.gene = self.result["gene_summary"]["input_genes"][0]

    def test_required_fields(self):
        expected = {"rank", "gene", "combined_D", "prior", "log_bf",
                    "combined", "N", "location", "status"}
        self.assertTrue(expected.issubset(set(self.gene.keys())))

    def test_no_unverified_fields(self):
        """Gene entries should NOT contain unverified columns."""
        unverified = {"prior_adj", "combined_adj", "positive_control"}
        self.assertFalse(unverified & set(self.gene.keys()))

    def test_combined_D_in_range(self):
        self.assertGreaterEqual(self.gene["combined_D"], 0)
        self.assertLessEqual(self.gene["combined_D"], 1)

    def test_location_format(self):
        loc = self.gene["location"]
        # Should be "Chrom:Start-End" format
        self.assertIn(":", loc)
        self.assertIn("-", loc)


class TestGeneSetSummary(unittest.TestCase):
    """Test gene-set summary ranking and fields."""

    def setUp(self):
        self.result = _get_full_interpretation()
        self.gss = self.result["gene_set_summary"]

    def test_top_gene_sets_present(self):
        self.assertTrue(len(self.gss["top_gene_sets"]) > 0)

    def test_sorted_by_abs_beta(self):
        betas = [abs(gs["beta"]) for gs in self.gss["top_gene_sets"]]
        self.assertEqual(betas, sorted(betas, reverse=True))

    def test_gene_set_fields(self):
        gs = self.gss["top_gene_sets"][0]
        expected = {"rank", "gene_set", "label", "N", "beta", "avg_postp", "P"}
        self.assertTrue(expected.issubset(set(gs.keys())))

    def test_total_gene_sets(self):
        self.assertEqual(self.gss["total_gene_sets"], 5)


class TestGenePathwayLinks(unittest.TestCase):
    """Test gene-pathway links from ggss.out."""

    def setUp(self):
        self.result = _get_full_interpretation()
        self.links = self.result["gene_pathway_links"]

    def test_by_gene_populated(self):
        self.assertTrue(len(self.links["by_gene"]) > 0)

    def test_gene_a_has_links(self):
        gene_a = self.links["by_gene"].get("GENE_A")
        self.assertIsNotNone(gene_a)
        self.assertTrue(len(gene_a["contributing_gene_sets"]) > 0)

    def test_gene_a_status(self):
        gene_a = self.links["by_gene"]["GENE_A"]
        self.assertEqual(gene_a["status"], "INPUT")

    def test_candidate_gene_links(self):
        gene_d = self.links["by_gene"].get("GENE_D")
        self.assertIsNotNone(gene_d)
        self.assertEqual(gene_d["status"], "CANDIDATE")

    def test_contributing_gene_sets_have_beta(self):
        gene_a = self.links["by_gene"]["GENE_A"]
        for gs_link in gene_a["contributing_gene_sets"]:
            self.assertIn("beta", gs_link)
            self.assertIn("gene_set", gs_link)


class TestInterpretationContracts(unittest.TestCase):
    """Test interpretation contracts are present and complete."""

    def test_all_verified_metrics_have_contracts(self):
        expected_keys = {"combined_D", "prior", "log_bf", "combined",
                         "beta_gss", "avg_postp", "P_gss"}
        self.assertEqual(set(INTERPRETATION_CONTRACTS.keys()), expected_keys)

    def test_contract_fields(self):
        for name, contract in INTERPRETATION_CONTRACTS.items():
            self.assertIn("definition", contract, f"Missing definition for {name}")
            self.assertIn("scale", contract, f"Missing scale for {name}")
            self.assertIn("anti_interpretation", contract,
                          f"Missing anti_interpretation for {name}")
            self.assertIn("stochastic", contract, f"Missing stochastic for {name}")

    def test_combined_D_description_model_derived(self):
        """combined_D must always be described as 'model-derived'."""
        defn = INTERPRETATION_CONTRACTS["combined_D"]["definition"]
        self.assertIn("model-derived", defn.lower())

    def test_combined_D_anti_interpretation(self):
        anti = INTERPRETATION_CONTRACTS["combined_D"]["anti_interpretation"]
        self.assertIn("NOT", anti)


class TestForbiddenLanguageScan(unittest.TestCase):
    """Scan ALL interpretation output strings for forbidden terms."""

    def setUp(self):
        self.result = _get_full_interpretation()

    def _collect_strings(self, obj, path=""):
        """Recursively collect all string values from a nested structure."""
        strings = []
        if isinstance(obj, str):
            strings.append((path, obj))
        elif isinstance(obj, dict):
            for k, v in obj.items():
                strings.extend(self._collect_strings(v, f"{path}.{k}"))
        elif isinstance(obj, (list, tuple)):
            for i, v in enumerate(obj):
                strings.extend(self._collect_strings(v, f"{path}[{i}]"))
        return strings

    def test_no_forbidden_terms_in_output(self):
        all_strings = self._collect_strings(self.result)
        for path, text in all_strings:
            text_lower = text.lower()
            for term in FORBIDDEN_TERMS:
                # Allow "anti_interpretation" fields to mention forbidden terms
                # as part of what NOT to say
                if "anti_interpretation" in path:
                    continue
                self.assertNotIn(
                    term.lower(), text_lower,
                    f"Forbidden term '{term}' found at {path}: '{text[:80]}'"
                )

    def test_forbidden_terms_list_nonempty(self):
        self.assertTrue(len(FORBIDDEN_TERMS) >= 7)


class TestConvergenceNote(unittest.TestCase):
    """Test stability note generation (post-audit terminology)."""

    def test_stability_met_note(self):
        result = _get_full_interpretation()
        note = result["convergence_note"]
        self.assertIn("PIGEAN stability criterion", note)
        self.assertIn("MET", note)
        self.assertIn("suitable for interpretation", note)
        self.assertIn("Formal MCMC convergence was not assessed", note)

    def test_none_convergence(self):
        gs = parse_gene_stats(os.path.join(FIXTURES, "gs_positive_controls.out"))
        result = interpret_results(gs, None, None, None, None,
                                   "positive-controls")
        self.assertIn("not assessed", result["convergence_note"].lower())

    def test_no_bare_converged_in_note(self):
        """The convergence note must NOT say bare 'CONVERGED'."""
        result = _get_full_interpretation()
        note = result["convergence_note"]
        # Should NOT contain bare "CONVERGED" without "stability criterion"
        self.assertNotIn("Model convergence: CONVERGED", note)


class TestCaveats(unittest.TestCase):
    """Test caveat generation."""

    def setUp(self):
        self.result = _get_full_interpretation()

    def test_stochastic_caveat_present(self):
        caveats_text = " ".join(self.result["caveats"])
        self.assertIn("stochastic", caveats_text.lower())

    def test_positive_control_caveat_present(self):
        """Positive-controls analysis must warn about input genes."""
        caveats_text = " ".join(self.result["caveats"])
        self.assertIn("input genes", caveats_text.lower())

    def test_model_not_causal_caveat(self):
        """Must mention that high scores mean model-prioritized, not causal."""
        caveats_text = " ".join(self.result["caveats"])
        self.assertIn("model", caveats_text.lower())


class TestNoneInputs(unittest.TestCase):
    """Test graceful handling of None inputs."""

    def test_all_none(self):
        result = interpret_results(None, None, None, None, None,
                                   "positive-controls")
        self.assertIsNone(result["gene_summary"])
        self.assertIsNone(result["gene_set_summary"])
        self.assertIsNone(result["gene_pathway_links"])

    def test_gs_only(self):
        gs = parse_gene_stats(os.path.join(FIXTURES, "gs_positive_controls.out"))
        result = interpret_results(gs, None, None, None, None,
                                   "positive-controls")
        self.assertIsNotNone(result["gene_summary"])
        self.assertIsNone(result["gene_set_summary"])
        self.assertIsNone(result["gene_pathway_links"])


if __name__ == "__main__":
    unittest.main()
