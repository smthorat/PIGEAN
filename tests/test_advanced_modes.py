"""Phase 6 tests for high-level model-mode selection."""

import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.config import resolve_config
from pigean.engine import build_priors_command
from pigean.modes import (
    COMPATIBILITY_MATRIX,
    MODES,
    SUPPORTED_MODES,
    get_mode,
)
from pigean.modes.base import CompatibilityStatus


class TestModeRegistry(unittest.TestCase):
    def test_only_verified_modes_registered(self):
        self.assertEqual(set(SUPPORTED_MODES), {"standard", "naive-priors"})
        self.assertNotIn("factor", MODES)
        self.assertNotIn("phewas", MODES)

    def test_unknown_mode_never_falls_back(self):
        with self.assertRaises(ValueError):
            get_mode("factor")

    def test_default_is_standard(self):
        config = resolve_config({})
        self.assertEqual(config["mode"], "standard")

    def test_explicit_mode_resolves(self):
        config = resolve_config({"mode": "naive-priors"})
        self.assertEqual(config["mode"], "naive-priors")

    def test_mode_metadata_is_manifest_ready(self):
        metadata = get_mode("naive-priors").describe({
            "analysis": "positive-controls",
        })
        self.assertEqual(metadata["name"], "naive-priors")
        self.assertEqual(metadata["engine_arguments"], ["naive_priors"])
        self.assertEqual(metadata["compatibility_status"], "SUPPORTED")
        self.assertEqual(metadata["stability_applicability"], "NOT_APPLICABLE")
        self.assertEqual(
            metadata["interpretation_compatibility"],
            "PARTIALLY_COMPATIBLE",
        )


class TestCompatibilityMatrix(unittest.TestCase):
    def test_standard_positive_controls_supported(self):
        self.assertEqual(
            COMPATIBILITY_MATRIX["positive-controls"]["standard"],
            "SUPPORTED",
        )

    def test_naive_positive_controls_supported(self):
        self.assertEqual(
            COMPATIBILITY_MATRIX["positive-controls"]["naive-priors"],
            "SUPPORTED",
        )

    def test_existing_broken_modes_are_engine_blocked(self):
        for evidence in ("gene-z-score", "gene-percentile"):
            self.assertEqual(
                COMPATIBILITY_MATRIX[evidence]["standard"],
                "ENGINE_BLOCKED",
            )
            self.assertEqual(
                COMPATIBILITY_MATRIX[evidence]["naive-priors"],
                "ENGINE_BLOCKED",
            )

    def test_unverified_advanced_combinations_rejected(self):
        mode = get_mode("naive-priors")
        for evidence in ("gene-bayes-factor", "exome", "gwas"):
            self.assertEqual(
                mode.compatibility_for(evidence),
                CompatibilityStatus.NOT_VERIFIED,
            )
            self.assertTrue(mode.validate_config({"analysis": evidence}))

    def test_internal_advanced_parameters_are_rejected(self):
        mode = get_mode("standard")
        for parameter in ("anchor", "phi", "alpha0"):
            issues = mode.validate_config({
                "analysis": "positive-controls",
                parameter: "supplied",
            })
            self.assertTrue(any(parameter in issue for issue in issues))


class TestEngineSubcommand(unittest.TestCase):
    def _build(self, subcommand):
        return build_priors_command(
            config={
                "max_num_gene_sets": 5000,
                "debug_level": 3,
                "gene_filter_value": 1,
                "gene_set_filter_value": 0.01,
            },
            evidence_engine_args=["--positive-controls-in", "/input.txt"],
            reference_paths={
                "gene_loc": "/gene.loc",
                "tss_loc": "/tss.loc",
                "exons_loc": "/exons.loc",
            },
            gene_set_paths=["/sets.txt"],
            gene_map_path="/map.txt",
            output_dir="/out",
            base_dir="/app",
            engine_subcommand=subcommand,
        )

    def test_standard_uses_gibbs(self):
        self.assertEqual(self._build("gibbs")[3], "gibbs")

    def test_naive_uses_exact_engine_token(self):
        self.assertEqual(self._build("naive_priors")[3], "naive_priors")


class TestInvalidModeCli(unittest.TestCase):
    def test_unknown_mode_writes_fail_manifest_without_engine(self):
        repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        with tempfile.TemporaryDirectory() as tmp:
            output = os.path.join(tmp, "run")
            proc = subprocess.run(
                [
                    sys.executable,
                    os.path.join(repo, "run_pigean.py"),
                    "--mode", "factor",
                    "--input", os.path.join(repo, "examples", "gene_list"),
                    "--output", output,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            with open(os.path.join(output, "run_manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertEqual(manifest["validation"]["status"], "FAIL")
            self.assertEqual(manifest["execution"]["status"], "NOT_RUN")
            self.assertEqual(manifest["advanced_mode"]["name"], "factor")


if __name__ == "__main__":
    unittest.main()
