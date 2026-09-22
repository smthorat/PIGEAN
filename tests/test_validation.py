"""Tests for pigean/validation.py and pigean/adapters/positive_controls.py."""
import sys, os, tempfile, shutil, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.adapters import ValidationStatus
from pigean.adapters.positive_controls import PositiveControlsAdapter
from pigean.validation import (
    load_gene_map, validate_input_genes, compute_validation_status
)
from pigean.adapters import (
    ValidationStatus, FileValidationResult, GeneQCResult,
    GeneSetPreparationResult, BackgroundQCResult,
)

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
GENE_MAP = os.path.join(REPO_ROOT, "data", "portal_gencode.gene.map")
GENE_LOC = os.path.join(REPO_ROOT, "data", "NCBI37.3.plink.gene.loc")
GENE_SETS = [
    os.path.join(REPO_ROOT, "data", "gene_set_list_mouse_2024.txt"),
    os.path.join(REPO_ROOT, "data", "gene_set_list_msigdb_nohp.txt"),
]

class TestFileValidation(unittest.TestCase):
    def setUp(self):
        self.adapter = PositiveControlsAdapter()
        self.tmpdir = tempfile.mkdtemp()
        self.params = {}

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_valid_input(self):
        path = os.path.join(self.tmpdir, "genes.txt")
        with open(path, "w") as f:
            f.write("LEP\nADIPOQ\nBRD2\n")
        result = self.adapter.validate_file(path, self.params)
        self.assertEqual(result.status, ValidationStatus.PASS)

    def test_nonexistent_input(self):
        result = self.adapter.validate_file("/does/not/exist", self.params)
        self.assertEqual(result.status, ValidationStatus.FAIL)
        self.assertTrue(any("not found" in i.lower() or "does not exist" in i.lower() or "not exist" in i.lower() for i in result.issues))

    def test_empty_input(self):
        path = os.path.join(self.tmpdir, "empty.txt")
        with open(path, "w") as f:
            pass
        result = self.adapter.validate_file(path, self.params)
        self.assertEqual(result.status, ValidationStatus.FAIL)

    def test_whitespace_only_input(self):
        path = os.path.join(self.tmpdir, "whitespace.txt")
        with open(path, "w") as f:
            f.write("   \n\n  \n")
        result = self.adapter.validate_file(path, self.params)
        # Should fail or at least not crash
        self.assertIn(result.status, [ValidationStatus.FAIL, ValidationStatus.PASS_WITH_WARNINGS])

class TestNormalization(unittest.TestCase):
    def setUp(self):
        self.adapter = PositiveControlsAdapter()
        self.tmpdir = tempfile.mkdtemp()
        self.params = {}

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_normal_input(self):
        path = os.path.join(self.tmpdir, "genes.txt")
        with open(path, "w") as f:
            f.write("LEP\nADIPOQ\nBRD2\n")
        outdir = os.path.join(self.tmpdir, "output")
        os.makedirs(outdir)
        result = self.adapter.normalize(path, outdir, self.params)
        self.assertEqual(result.original_count, 3)
        self.assertEqual(result.normalized_count, 3)
        self.assertTrue(os.path.exists(os.path.join(outdir, "input", "positive_controls.txt")))
        self.assertTrue(os.path.exists(result.normalized_path))

    def test_duplicates_removed(self):
        path = os.path.join(self.tmpdir, "genes.txt")
        with open(path, "w") as f:
            f.write("LEP\nLEP\nADIPOQ\n")
        outdir = os.path.join(self.tmpdir, "output")
        os.makedirs(outdir)
        result = self.adapter.normalize(path, outdir, self.params)
        self.assertEqual(result.duplicates_removed, 1)
        self.assertEqual(result.normalized_count, 2)

    def test_blank_lines_removed(self):
        path = os.path.join(self.tmpdir, "genes.txt")
        with open(path, "w") as f:
            f.write("LEP\n\n\nADIPOQ\n  \nBRD2\n")
        outdir = os.path.join(self.tmpdir, "output")
        os.makedirs(outdir)
        result = self.adapter.normalize(path, outdir, self.params)
        self.assertEqual(result.normalized_count, 3)
        self.assertGreater(result.blanks_removed, 0)

    def test_whitespace_stripped(self):
        path = os.path.join(self.tmpdir, "genes.txt")
        with open(path, "w") as f:
            f.write("  LEP  \n  ADIPOQ\n")
        outdir = os.path.join(self.tmpdir, "output")
        os.makedirs(outdir)
        result = self.adapter.normalize(path, outdir, self.params)
        # Read normalized output
        with open(result.normalized_path) as f:
            lines = [l.strip() for l in f if l.strip()]
        self.assertEqual(lines, ["LEP", "ADIPOQ"])

class TestGeneQC(unittest.TestCase):
    """These tests require the actual reference files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Only run if reference files exist
        if not os.path.exists(GENE_MAP):
            self.skipTest("Reference files not available")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_known_genes_recognized(self):
        gene_map = load_gene_map(GENE_MAP)
        self.assertIn("LEP", gene_map)
        self.assertIn("ADIPOQ", gene_map)
        self.assertIn("BRD2", gene_map)

    def test_fake_genes_unrecognized(self):
        gene_map = load_gene_map(GENE_MAP)
        self.assertNotIn("FAKEGENE1", gene_map)
        self.assertNotIn("NOTREAL", gene_map)

    def test_all_fake_genes_fail(self):
        """Zero recognized genes → FAIL."""
        norm_path = os.path.join(self.tmpdir, "fake.txt")
        with open(norm_path, "w") as f:
            f.write("FAKEGENE1\n")
        qc = validate_input_genes(norm_path, GENE_MAP, GENE_LOC, GENE_SETS)
        self.assertEqual(qc.recognized_count, 0)
        self.assertEqual(qc.status, ValidationStatus.FAIL)

    def test_mixed_valid_invalid(self):
        norm_path = os.path.join(self.tmpdir, "mixed.txt")
        with open(norm_path, "w") as f:
            f.write("LEP\nFAKEGENE1\nADIPOQ\n")
        qc = validate_input_genes(norm_path, GENE_MAP, GENE_LOC, GENE_SETS)
        self.assertEqual(qc.recognized_count, 2)
        self.assertEqual(qc.unresolved_count, 1)
        self.assertEqual(qc.status, ValidationStatus.PASS_WITH_WARNINGS)

    def test_plink1_unresolved(self):
        """PLINK1 from the standard input should be unresolved."""
        norm_path = os.path.join(self.tmpdir, "genes.txt")
        with open(norm_path, "w") as f:
            f.write("LEP\nPLINK1\nADIPOQ\n")
        qc = validate_input_genes(norm_path, GENE_MAP, GENE_LOC, GENE_SETS)
        # Find PLINK1 in per_gene
        plink1 = [g for g in qc.per_gene if g["gene"] == "PLINK1"]
        self.assertEqual(len(plink1), 1)
        self.assertFalse(plink1[0]["recognized"])

class TestCombinedStatus(unittest.TestCase):
    def test_both_pass(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        self.assertEqual(compute_validation_status(fr, qr), ValidationStatus.PASS)

    def test_file_fail(self):
        fr = FileValidationResult(status=ValidationStatus.FAIL, issues=["bad file"])
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        self.assertEqual(compute_validation_status(fr, qr), ValidationStatus.FAIL)

    def test_qc_warn(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS_WITH_WARNINGS, input_count=5, recognized_count=3, unresolved_count=2)
        self.assertEqual(compute_validation_status(fr, qr), ValidationStatus.PASS_WITH_WARNINGS)


class TestCombinedStatusPhase2(unittest.TestCase):
    """Phase 2 tests for compute_validation_status with new parameters."""

    def test_background_file_fail_causes_fail(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        bg_fr = FileValidationResult(status=ValidationStatus.FAIL, issues=["bad background"])
        self.assertEqual(
            compute_validation_status(fr, qr, background_file_result=bg_fr),
            ValidationStatus.FAIL
        )

    def test_background_qc_warning_propagates(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        bg_qc = BackgroundQCResult(
            status=ValidationStatus.PASS_WITH_WARNINGS,
            issues=["1 input gene not in background"],
            total_count=10, mapped_count=10, unmapped_count=0,
            input_genes_missing_from_background=["MISSING1"],
        )
        self.assertEqual(
            compute_validation_status(fr, qr, background_qc_result=bg_qc),
            ValidationStatus.PASS_WITH_WARNINGS
        )

    def test_background_qc_fail_causes_fail(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        bg_qc = BackgroundQCResult(
            status=ValidationStatus.FAIL,
            issues=["No recognized genes"],
            total_count=2, mapped_count=0, unmapped_count=2,
        )
        self.assertEqual(
            compute_validation_status(fr, qr, background_qc_result=bg_qc),
            ValidationStatus.FAIL
        )

    def test_custom_gs_fail_causes_fail(self):
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        gs_prep = GeneSetPreparationResult(
            status=ValidationStatus.FAIL,
            issues=["missing files"],
        )
        self.assertEqual(
            compute_validation_status(fr, qr, custom_gs_prep_result=gs_prep),
            ValidationStatus.FAIL
        )

    def test_all_none_optional_same_as_phase1(self):
        """All optional params None → same behavior as Phase 1."""
        fr = FileValidationResult(status=ValidationStatus.PASS)
        qr = GeneQCResult(status=ValidationStatus.PASS, input_count=5, recognized_count=5)
        self.assertEqual(
            compute_validation_status(fr, qr, None, None, None),
            ValidationStatus.PASS
        )

if __name__ == "__main__":
    unittest.main()
