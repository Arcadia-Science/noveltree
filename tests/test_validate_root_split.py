import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "validate_root_split.py"
SPEC = importlib.util.spec_from_file_location("validate_root_split", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ValidateRootSplitTests(unittest.TestCase):
    def test_accepts_child_order_change_and_records_final_validation(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            expected = root / "expected.nwk"
            observed = root / "observed.nwk"
            qc = root / "qc.tsv"
            expected.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
            observed.write_text("((D:2,C:2):2,(B:2,A:2):2);\n")
            qc.write_text("metric\tvalue\nroot_split_validated\ttrue\n")

            MODULE.validate(expected, observed, qc)

            self.assertIn("final_root_split_validated\ttrue", qc.read_text())

    def test_rejects_a_changed_root(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            expected = root / "expected.nwk"
            observed = root / "observed.nwk"
            expected.write_text("((A,B),(C,D));\n")
            observed.write_text("((A,C),(B,D));\n")
            with self.assertRaisesRegex(ValueError, "root split differs"):
                MODULE.validate(expected, observed, root / "qc.tsv")


if __name__ == "__main__":
    unittest.main()
