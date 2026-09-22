"""Tests for pigean/config.py — config resolution and precedence."""
import sys, os, json, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pigean.config import DEFAULTS, resolve_config

class MockArgs:
    """Mock argparse Namespace with all None values."""
    def __init__(self, **kwargs):
        # Start with all None
        for key in ["analysis", "gene_sets", "genome_build", "preset",
                     "max_num_gene_sets", "gene_filter_value",
                     "gene_set_filter_value", "debug_level",
                     "custom_gene_set_format", "custom_gene_set_action"]:
            setattr(self, key, None)
        # Override with provided values
        for k, v in kwargs.items():
            setattr(self, k, v)

class TestDefaults(unittest.TestCase):
    def test_defaults_match_phase0(self):
        """DEFAULTS must match the Phase 0 scientific parameters."""
        self.assertEqual(DEFAULTS["analysis"], "positive-controls")
        self.assertEqual(DEFAULTS["genome_build"], "hg19")
        self.assertEqual(DEFAULTS["gene_sets"], "default")
        self.assertEqual(DEFAULTS["max_num_gene_sets"], 5000)
        self.assertEqual(DEFAULTS["gene_filter_value"], 1)
        self.assertEqual(DEFAULTS["gene_set_filter_value"], 0.01)
        self.assertEqual(DEFAULTS["debug_level"], 3)

class TestConfigPrecedence(unittest.TestCase):
    def test_defaults_only(self):
        args = MockArgs()
        config = resolve_config(args)
        self.assertEqual(config["max_num_gene_sets"], 5000)
        self.assertEqual(config["gene_filter_value"], 1)

    def test_json_overrides_defaults(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"max_num_gene_sets": 3000}, f)
            f.flush()
            args = MockArgs()
            config = resolve_config(args, config_file=f.name)
            self.assertEqual(config["max_num_gene_sets"], 3000)
            # Other defaults preserved
            self.assertEqual(config["debug_level"], 3)
        os.unlink(f.name)

    def test_cli_overrides_json(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"max_num_gene_sets": 3000}, f)
            f.flush()
            args = MockArgs(max_num_gene_sets=2000)
            config = resolve_config(args, config_file=f.name)
            self.assertEqual(config["max_num_gene_sets"], 2000)
        os.unlink(f.name)

    def test_none_cli_does_not_override(self):
        """CLI args that are None must not overwrite config/defaults."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"gene_sets": "default"}, f)
            f.flush()
            args = MockArgs()  # gene_sets=None
            config = resolve_config(args, config_file=f.name)
            self.assertEqual(config["gene_sets"], "default")
        os.unlink(f.name)


class TestPhase2Config(unittest.TestCase):
    """Tests for Phase 2 config keys (custom gene sets)."""

    def test_custom_gene_set_format_from_cli(self):
        """custom_gene_set_format from CLI is picked up."""
        args = MockArgs(custom_gene_set_format="gmt")
        config = resolve_config(args)
        self.assertEqual(config.get("custom_gene_set_format"), "gmt")

    def test_custom_gene_set_action_from_cli(self):
        """custom_gene_set_action from CLI is picked up."""
        args = MockArgs(custom_gene_set_action="replace")
        config = resolve_config(args)
        self.assertEqual(config.get("custom_gene_set_action"), "replace")

    def test_custom_gene_set_action_from_json(self):
        """custom_gene_set_action from JSON config works."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"custom_gene_set_action": "replace"}, f)
            f.flush()
            args = MockArgs()
            config = resolve_config(args, config_file=f.name)
            self.assertEqual(config.get("custom_gene_set_action"), "replace")
        os.unlink(f.name)

    def test_cli_overrides_json_for_custom_action(self):
        """CLI --custom-gene-set-action overrides JSON config."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"custom_gene_set_action": "replace"}, f)
            f.flush()
            args = MockArgs(custom_gene_set_action="supplement")
            config = resolve_config(args, config_file=f.name)
            self.assertEqual(config["custom_gene_set_action"], "supplement")
        os.unlink(f.name)

    def test_none_custom_keys_not_in_defaults(self):
        """Phase 2 keys with None CLI values don't appear in resolved config
        unless set by JSON config."""
        args = MockArgs()
        config = resolve_config(args)
        # These keys are not in DEFAULTS, so with None CLI they're absent
        self.assertNotIn("custom_gene_set_format", config)
        self.assertNotIn("custom_gene_set_action", config)

if __name__ == "__main__":
    unittest.main()
