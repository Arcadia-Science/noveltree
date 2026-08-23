import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "resolve_polytomies.py"


class ResolvePolytomiesTests(unittest.TestCase):
    def run_script(self, root, tree_text, *extra_args):
        source = root / "input.newick"
        output = root / "output.newick"
        source.write_text(tree_text)
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                str(source),
                str(output),
                *extra_args,
            ],
            capture_output=True,
            text=True,
        )
        return result, output

    def test_default_behavior_keeps_zero_length_resolution_edges(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result, output = self.run_script(root, "(A:1,B:2,C:3)root;\n")

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_text(), "(A:1,(B:2,C:3):0)root;\n")

    def test_clamps_every_non_root_branch_after_resolution(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            qc = root / "qc.tsv"
            result, output = self.run_script(
                root,
                "((A:0,B:9e-9,C:-2):nan,D:101,E)root:0;\n",
                "--min-branch-length",
                "1e-6",
                "--max-branch-length",
                "100",
                "--branch-length-qc",
                str(qc),
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                output.read_text(),
                "((A:1e-06,(B:1e-06,C:1e-06):1e-06):1e-06,"
                "(D:100,E:1e-06):1e-06)root:0;\n",
            )
            self.assertIn("8 checked; 8 adjusted", result.stdout)

            with qc.open() as handle:
                rows = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["branches_total"], "8")
            self.assertEqual(rows["branches_adjusted"], "8")
            self.assertEqual(rows["branches_missing"], "1")
            self.assertEqual(rows["branches_nonfinite_or_invalid"], "1")
            self.assertEqual(rows["branches_negative"], "1")
            # One input zero plus two zero-length edges introduced while
            # resolving the two three-way polytomies.
            self.assertEqual(rows["branches_zero"], "3")
            self.assertEqual(rows["branches_positive_below_minimum"], "1")
            self.assertEqual(rows["branches_above_maximum"], "1")
            self.assertEqual(rows["minimum_allowed"], "1e-06")
            self.assertEqual(rows["maximum_allowed"], "100")

    def test_rejects_incomplete_or_invalid_bounds(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result, _output = self.run_script(
                root,
                "(A:1,B:1);\n",
                "--min-branch-length",
                "1e-6",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("must be supplied together", result.stderr)

            result, _output = self.run_script(
                root,
                "(A:1,B:1);\n",
                "--min-branch-length",
                "2",
                "--max-branch-length",
                "1",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("greater than or equal", result.stderr)

    def test_generax_enables_bounds_after_polytomy_resolution(self):
        module = (REPO / "modules" / "local" / "generax_per_species.nf").read_text()

        self.assertIn("--min-branch-length 1e-6", module)
        self.assertIn("--max-branch-length 100", module)
        self.assertIn("--branch-length-qc ${og}_branch_length_qc.tsv", module)
        self.assertIn("cp ${og}_branch_length_qc.tsv $og/", module)

        output_block = module.split("output:", 1)[1].split("when:", 1)[0]
        self.assertNotIn("branch_length_qc", output_block)


if __name__ == "__main__":
    unittest.main()
