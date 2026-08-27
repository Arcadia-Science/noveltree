import argparse
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))


def load_script(name):
    path = REPO / "bin" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PREPARE = load_script("prepare_reconciliation_tree.py")


class PrepareReconciliationTreeTests(unittest.TestCase):
    def test_resolves_polytomy_clamps_lengths_and_preserves_leaves(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = root / "tree.nwk"
            tree.write_text(
                "(Species-a_p1:0,Species-b_p2:1e-9,Species-c_p3:101);\n"
            )
            alignment = root / "alignment.fa"
            alignment.write_text(
                ">Species-a_p1\nAA\n>Species-b_p2\nAA\n>Species-c_p3\nAA\n"
            )
            mapping = root / "map.link"
            mapping.write_text(
                "Species-a_p1\tSpecies-a\n"
                "Species-b_p2\tSpecies-b\n"
                "Species-c_p3\tSpecies-c\n"
            )
            output = root / "ready.nwk"
            qc = root / "qc.tsv"
            PREPARE.prepare(
                argparse.Namespace(
                    orthogroup="OG1",
                    tree=tree,
                    alignment=alignment,
                    mapping=mapping,
                    output=output,
                    qc=qc,
                    min_branch_length=1e-6,
                    max_branch_length=100.0,
                )
            )
            ready = PREPARE.parse_newick(output.read_text())
            self.assertEqual(
                set(PREPARE.leaf_labels(ready)),
                {"Species-a_p1", "Species-b_p2", "Species-c_p3"},
            )
            with qc.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["multifurcations_resolved"], "1")
            self.assertEqual(row["branches_adjusted"], "4")
            self.assertEqual(row["branches_zero"], "2")
            self.assertEqual(row["branches_positive_below_minimum"], "1")
            self.assertEqual(row["branches_above_maximum"], "1")

    def test_rejects_alignment_or_mapping_leaf_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = root / "tree.nwk"
            tree.write_text("(Species-a_p1:1,Species-b_p2:1);\n")
            alignment = root / "alignment.fa"
            alignment.write_text(">Species-a_p1\nAA\n")
            mapping = root / "map.link"
            mapping.write_text(
                "Species-a_p1\tSpecies-a\nSpecies-b_p2\tSpecies-b\n"
            )
            with self.assertRaisesRegex(ValueError, "tree/alignment leaf mismatch"):
                PREPARE.prepare(
                    argparse.Namespace(
                        orthogroup="OG1",
                        tree=tree,
                        alignment=alignment,
                        mapping=mapping,
                        output=root / "ready.nwk",
                        qc=root / "qc.tsv",
                        min_branch_length=1e-6,
                        max_branch_length=100.0,
                    )
                )


if __name__ == "__main__":
    unittest.main()
