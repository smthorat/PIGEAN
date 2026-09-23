"""Phase 6 tests for the verified naive-priors wrapper mode."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.convergence import (
    FormalConvergenceStatus,
    IterationStatus,
    StabilityStatus,
    make_not_applicable_result,
    write_convergence_json,
)
from pigean.modes import get_mode
from pigean.modes.base import (
    InterpretationCompatibility,
    StabilityApplicability,
)
from pigean.parsers import parse_gene_stats
from pigean.interpretation import interpret_results


class TestNaivePriorsMode(unittest.TestCase):
    def setUp(self):
        self.mode = get_mode("naive-priors")

    def test_engine_subcommand(self):
        self.assertEqual(self.mode.ENGINE_SUBCOMMAND, "naive_priors")

    def test_no_raw_advanced_flags(self):
        self.assertEqual(
            self.mode.build_engine_args({"analysis": "positive-controls"}),
            [],
        )

    def test_positive_controls_valid(self):
        self.assertEqual(
            self.mode.validate_config({"analysis": "positive-controls"}), []
        )

    def test_trace_rejected(self):
        issues = self.mode.validate_config({
            "analysis": "positive-controls",
            "enable_convergence_trace": True,
        })
        self.assertTrue(any("incompatible" in issue for issue in issues))

    def test_stability_not_applicable(self):
        self.assertEqual(
            self.mode.STABILITY_APPLICABILITY,
            StabilityApplicability.NOT_APPLICABLE,
        )

    def test_interpretation_partially_compatible(self):
        self.assertEqual(
            self.mode.INTERPRETATION_COMPATIBILITY,
            InterpretationCompatibility.PARTIALLY_COMPATIBLE,
        )


class TestNaiveStabilityResult(unittest.TestCase):
    def test_not_applicable_statuses(self):
        result = make_not_applicable_result("No outer Gibbs loop.")
        self.assertEqual(
            result["stability_status"], StabilityStatus.NOT_APPLICABLE
        )
        self.assertEqual(
            result["formal_convergence_status"],
            FormalConvergenceStatus.NOT_APPLICABLE,
        )
        self.assertEqual(
            result["iteration_status"], IterationStatus.NOT_APPLICABLE
        )

    def test_json_serialization(self):
        result = make_not_applicable_result("No outer Gibbs loop.")
        with tempfile.TemporaryDirectory() as tmp:
            path = write_convergence_json(result, tmp)
            with open(path) as fh:
                data = json.load(fh)
        self.assertEqual(data["stability_status"], "NOT_APPLICABLE")
        self.assertEqual(data["formal_convergence_status"], "NOT_APPLICABLE")


class TestNaiveOutputParsing(unittest.TestCase):
    def test_combined_fallback_ranking_without_combined_d(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "gs.out")
            with open(path, "w") as fh:
                fh.write("Gene\tprior\tcombined\tlog_bf\tN\n")
                fh.write("LOW\t0.1\t0.5\t0.4\t2\n")
                fh.write("HIGH\t0.2\t1.5\t1.3\t3\n")
            parsed = parse_gene_stats(path)
        self.assertEqual(parsed["genes"][0]["Gene"], "HIGH")

    def test_interpretation_omits_combined_d_contract(self):
        parsed = {
            "genes": [{
                "Gene": "GENE1", "prior": "0.2", "combined": "1.2",
                "log_bf": "1.0", "N": "3",
            }],
            "row_count": 1,
            "has_positive_control": False,
        }
        result = interpret_results(
            parsed, None, None, None,
            make_not_applicable_result("No outer Gibbs loop."),
            "positive-controls", advanced_mode="naive-priors",
        )
        self.assertEqual(result["ranking_metric"], "combined")
        self.assertNotIn("combined_D", result["interpretation_contracts"])
        caveats = " ".join(result["caveats"])
        self.assertNotIn("high combined_D", caveats)


if __name__ == "__main__":
    unittest.main()
