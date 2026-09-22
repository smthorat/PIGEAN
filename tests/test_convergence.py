"""Tests for pigean/convergence.py — Phase 5 stability assessment.

Post-audit: Uses corrected terminology.  The engine's max-fractional-SEM
criterion is reported as PIGEAN stability, NOT formal MCMC convergence.
"""

import json
import os
import tempfile
import unittest

from pigean.convergence import (
    ConvergenceStatus,  # backward-compat alias
    FormalConvergenceStatus,
    IterationStatus,
    StabilityStatus,
    assess_convergence,
    parse_convergence_from_log,
    write_convergence_json,
)
from pigean.parsers import parse_params

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "outputs")


class TestStabilityStatusEnum(unittest.TestCase):
    """Tests for the StabilityStatus enum."""

    def test_all_values(self):
        expected = {
            "PIGEAN_STABILITY_CRITERION_MET",
            "PIGEAN_STABILITY_CRITERION_UNCERTAIN",
            "PIGEAN_STABILITY_CRITERION_NOT_MET",
            "NOT_ASSESSED",
            "TRACE_NOT_AVAILABLE",
            "ENGINE_FAILED",
        }
        actual = {s.value for s in StabilityStatus}
        self.assertEqual(actual, expected)

    def test_value_access(self):
        self.assertEqual(
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET.value,
            "PIGEAN_STABILITY_CRITERION_MET",
        )

    def test_backward_compat_alias(self):
        """ConvergenceStatus is an alias for StabilityStatus."""
        self.assertIs(ConvergenceStatus, StabilityStatus)

    def test_no_converged_value(self):
        """The enum must NOT contain a bare 'CONVERGED' value."""
        values = {s.value for s in StabilityStatus}
        self.assertNotIn("CONVERGED", values)


class TestFormalConvergenceStatus(unittest.TestCase):
    """Formal MCMC convergence is always NOT_FORMALLY_ASSESSED."""

    def test_only_not_assessed(self):
        self.assertEqual(
            FormalConvergenceStatus.NOT_FORMALLY_ASSESSED.value,
            "NOT_FORMALLY_ASSESSED",
        )


class TestIterationStatus(unittest.TestCase):
    """Iteration status is separate from stability status."""

    def test_cap_reached(self):
        self.assertEqual(
            IterationStatus.MAX_ITERATION_CAP_REACHED.value,
            "MAX_ITERATION_CAP_REACHED",
        )

    def test_stopped_before_cap(self):
        self.assertEqual(
            IterationStatus.STOPPED_BEFORE_CAP.value,
            "STOPPED_BEFORE_CAP",
        )


class TestAssessStabilityCriterionMet(unittest.TestCase):
    """Tier 1: Early stop (num_gibbs_iter=200), no restarts → MET."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_converged.out"))
        self.result = assess_convergence(parsed)

    def test_stability_status_met(self):
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET,
        )

    def test_assessment_source_summary(self):
        self.assertEqual(self.result["assessment_source"], "SUMMARY_BASED")

    def test_formal_convergence_not_assessed(self):
        self.assertEqual(
            self.result["formal_convergence_status"],
            FormalConvergenceStatus.NOT_FORMALLY_ASSESSED,
        )

    def test_formal_convergence_reason(self):
        self.assertIn("chain traces", self.result["formal_convergence_reason"])

    def test_iteration_status_before_cap(self):
        self.assertEqual(
            self.result["iteration_status"],
            IterationStatus.STOPPED_BEFORE_CAP,
        )

    def test_evidence_present(self):
        self.assertTrue(len(self.result["evidence"]) > 0)

    def test_no_warnings(self):
        self.assertEqual(len(self.result["warnings"]), 0)

    def test_details(self):
        d = self.result["details"]
        self.assertEqual(d["num_gibbs_iter"], 200)
        self.assertEqual(d["num_gibbs_restarts"], 0)
        self.assertFalse(d["hit_iteration_cap"])

    def test_threshold_in_details(self):
        self.assertEqual(self.result["details"]["max_frac_sem_threshold"], 0.01)

    def test_stability_method(self):
        self.assertEqual(
            self.result["stability_method"], "PIGEAN_MAX_FRACTIONAL_SEM")

    def test_legacy_status_alias(self):
        """Legacy 'status' key still works."""
        self.assertEqual(
            self.result["status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET,
        )


class TestAssessStabilityUncertainHitCap(unittest.TestCase):
    """Tier 1: num_gibbs_iter=499 (hit cap), no restarts → UNCERTAIN."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_hit_cap.out"))
        self.result = assess_convergence(parsed)

    def test_stability_status_uncertain(self):
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_UNCERTAIN,
        )

    def test_assessment_source_summary(self):
        self.assertEqual(self.result["assessment_source"], "SUMMARY_BASED")

    def test_iteration_status_cap_reached(self):
        self.assertEqual(
            self.result["iteration_status"],
            IterationStatus.MAX_ITERATION_CAP_REACHED,
        )

    def test_has_warnings(self):
        self.assertTrue(len(self.result["warnings"]) > 0)

    def test_hit_cap_in_details(self):
        self.assertTrue(self.result["details"]["hit_iteration_cap"])


class TestAssessStabilityNotMetRestarted(unittest.TestCase):
    """Tier 1: Hit cap + 2 restarts → NOT_MET."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_restarted.out"))
        self.result = assess_convergence(parsed)

    def test_stability_status_not_met(self):
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_NOT_MET,
        )

    def test_details_restarts(self):
        self.assertEqual(self.result["details"]["num_gibbs_restarts"], 2)

    def test_hit_cap_true(self):
        self.assertTrue(self.result["details"]["hit_iteration_cap"])

    def test_iteration_status_cap_reached(self):
        self.assertEqual(
            self.result["iteration_status"],
            IterationStatus.MAX_ITERATION_CAP_REACHED,
        )


class TestAssessStabilityNoneParams(unittest.TestCase):
    """None input → NOT_ASSESSED."""

    def test_none_returns_not_assessed(self):
        result = assess_convergence(None)
        self.assertEqual(
            result["stability_status"], StabilityStatus.NOT_ASSESSED)
        self.assertEqual(result["assessment_source"], "INSUFFICIENT")

    def test_empty_params_returns_not_assessed(self):
        result = assess_convergence({"params": {}})
        self.assertEqual(
            result["stability_status"], StabilityStatus.NOT_ASSESSED)


class TestLogParsingConverged(unittest.TestCase):
    """Test log parsing with explicit precision message."""

    def setUp(self):
        self.log_path = os.path.join(FIXTURES, "log_converged.log")
        self.log_info = parse_convergence_from_log(self.log_path)

    def test_precision_achieved(self):
        self.assertTrue(self.log_info["precision_achieved"])

    def test_final_max_ratio(self):
        self.assertAlmostEqual(self.log_info["final_max_ratio"], 0.00606)

    def test_r_hat_entries(self):
        self.assertTrue(len(self.log_info["r_hat_entries"]) > 0)

    def test_final_gibbs_iteration(self):
        self.assertEqual(self.log_info["final_gibbs_iteration"], 200)


class TestLogParsingHitCap(unittest.TestCase):
    """Test log parsing with SEM below threshold but no precision message."""

    def setUp(self):
        self.log_path = os.path.join(FIXTURES, "log_hit_cap.log")
        self.log_info = parse_convergence_from_log(self.log_path)

    def test_no_precision_message(self):
        self.assertFalse(self.log_info["precision_achieved"])

    def test_final_max_ratio_below_threshold(self):
        self.assertLess(self.log_info["final_max_ratio"], 0.01)

    def test_final_gibbs_iteration(self):
        self.assertEqual(self.log_info["final_gibbs_iteration"], 500)


class TestLogParsingNoStability(unittest.TestCase):
    """Test log parsing with SEM ratio above threshold."""

    def setUp(self):
        self.log_path = os.path.join(FIXTURES, "log_no_convergence.log")
        self.log_info = parse_convergence_from_log(self.log_path)

    def test_final_max_ratio_above_threshold(self):
        self.assertGreater(self.log_info["final_max_ratio"], 0.01)

    def test_no_precision_message(self):
        self.assertFalse(self.log_info["precision_achieved"])


class TestLogParsingMissing(unittest.TestCase):
    """Missing log returns None."""

    def test_missing_log(self):
        result = parse_convergence_from_log("/tmp/nonexistent.log")
        self.assertIsNone(result)

    def test_none_path(self):
        result = parse_convergence_from_log(None)
        self.assertIsNone(result)


class TestTier2UpgradeHitCapWithLog(unittest.TestCase):
    """Tier 2: Hit cap (Tier 1=UNCERTAIN) + log SEM < 0.01 → MET (TLOG_BASED)."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_hit_cap.out"))
        log_path = os.path.join(FIXTURES, "log_hit_cap.log")
        self.result = assess_convergence(parsed, log_path=log_path)

    def test_upgraded_to_met(self):
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET,
        )

    def test_tlog_based(self):
        self.assertEqual(self.result["assessment_source"], "TLOG_BASED")

    def test_log_parsed(self):
        self.assertTrue(self.result["details"]["log_parsed"])

    def test_final_sem_ratio(self):
        self.assertAlmostEqual(
            self.result["details"]["final_max_sem_ratio"], 0.00606)

    def test_formal_convergence_still_not_assessed(self):
        """Formal MCMC convergence must NOT change with log data."""
        self.assertEqual(
            self.result["formal_convergence_status"],
            FormalConvergenceStatus.NOT_FORMALLY_ASSESSED,
        )

    def test_iteration_status_still_cap(self):
        """Iteration status must remain cap-reached regardless of SEM."""
        self.assertEqual(
            self.result["iteration_status"],
            IterationStatus.MAX_ITERATION_CAP_REACHED,
        )


class TestTier2PrecisionMessage(unittest.TestCase):
    """Tier 2: Explicit precision message → MET (TLOG_BASED)."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_converged.out"))
        log_path = os.path.join(FIXTURES, "log_converged.log")
        self.result = assess_convergence(parsed, log_path=log_path)

    def test_met(self):
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET,
        )

    def test_tlog_based(self):
        self.assertEqual(self.result["assessment_source"], "TLOG_BASED")

    def test_precision_message_in_details(self):
        self.assertTrue(self.result["details"]["precision_achieved_message"])


class TestTier2LogAboveThreshold(unittest.TestCase):
    """Tier 2: Log shows SEM ratio 0.076 (above threshold), does not upgrade."""

    def setUp(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_hit_cap.out"))
        log_path = os.path.join(FIXTURES, "log_no_convergence.log")
        self.result = assess_convergence(parsed, log_path=log_path)

    def test_still_uncertain(self):
        # p.out: hit cap, no restarts → UNCERTAIN
        # log: max_ratio 0.076 > 0.01 → does NOT upgrade
        self.assertEqual(
            self.result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_UNCERTAIN,
        )

    def test_log_parsed(self):
        self.assertTrue(self.result["details"]["log_parsed"])


class TestWriteConvergenceJson(unittest.TestCase):
    """Test convergence.json output with new schema."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        parsed = parse_params(os.path.join(FIXTURES, "p_converged.out"))
        self.result = assess_convergence(parsed)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_writes_file(self):
        path = write_convergence_json(self.result, self.tmpdir)
        self.assertTrue(os.path.isfile(path))

    def test_valid_json(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertIn("stability_status", data)
        self.assertIn("assessment_source", data)
        self.assertIn("details", data)

    def test_stability_status_is_string(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertIsInstance(data["stability_status"], str)
        self.assertEqual(
            data["stability_status"], "PIGEAN_STABILITY_CRITERION_MET")

    def test_no_converged_in_status(self):
        """The status field must NOT be bare 'CONVERGED'."""
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertNotEqual(data["stability_status"], "CONVERGED")

    def test_formal_convergence_present(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertEqual(
            data["formal_convergence_status"], "NOT_FORMALLY_ASSESSED")

    def test_iteration_status_present(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertIn("iteration_status", data)
        self.assertEqual(data["iteration_status"], "STOPPED_BEFORE_CAP")

    def test_stability_method_present(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertEqual(
            data["stability_method"], "PIGEAN_MAX_FRACTIONAL_SEM")

    def test_threshold_present(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertEqual(data["threshold"], 0.01)

    def test_formal_convergence_reason(self):
        path = write_convergence_json(self.result, self.tmpdir)
        with open(path) as f:
            data = json.load(f)
        self.assertIn("chain traces", data["formal_convergence_reason"])


class TestThresholdUnchanged(unittest.TestCase):
    """The 0.01 threshold must remain the engine's original value."""

    def test_default_threshold(self):
        parsed = parse_params(os.path.join(FIXTURES, "p_converged.out"))
        result = assess_convergence(parsed)
        self.assertEqual(result["details"]["max_frac_sem_threshold"], 0.01)

    def test_calculation_unchanged(self):
        """The same numerical ratio must produce the same status."""
        parsed = parse_params(os.path.join(FIXTURES, "p_hit_cap.out"))
        log_path = os.path.join(FIXTURES, "log_hit_cap.log")
        result = assess_convergence(parsed, log_path=log_path)
        # 0.00606 < 0.01 → MET
        self.assertAlmostEqual(
            result["details"]["final_max_sem_ratio"], 0.00606)
        self.assertEqual(
            result["stability_status"],
            StabilityStatus.PIGEAN_STABILITY_CRITERION_MET,
        )


if __name__ == "__main__":
    unittest.main()
