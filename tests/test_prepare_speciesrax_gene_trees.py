import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))

from resolve_polytomies import parse_newick  # noqa: E402


def walk(node):
    yield node
    for child in node["children"]:
        yield from walk(child)


class PrepareSpeciesRaxGeneTreesTests(unittest.TestCase):
    def run_preparation(self, root):
        return subprocess.run(
            [
                sys.executable,
                str(REPO / "bin" / "prepare_speciesrax_gene_trees.py"),
                "--manifest-glob",
                "speciesrax_inputs_*.tsv",
                "--report",
                "speciesrax_gene_tree_validation.tsv",
            ],
            cwd=root,
            capture_output=True,
            text=True,
        )

    def write_manifest(self, root, families):
        manifest = root / "speciesrax_inputs_00000000.tsv"
        with manifest.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["orthogroup", "gene_tree", "mapping"])
            for orthogroup, tree in families:
                mapping = f"{orthogroup}_map.link"
                (root / mapping).write_text("gene species\n")
                writer.writerow([orthogroup, tree, mapping])

    def test_resolves_root_and_internal_polytomies_without_rewriting_binary_tree(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            binary_name = "OG0000001_ft.newick"
            polytomy_name = "OG0000002_ft.newick"
            binary = "(A:1,B:2)0.9:3;\n"
            (root / binary_name).write_text(binary)
            (root / polytomy_name).write_text(
                "((A:1,B:2,C:3)0.8:4,D:5,E:6)0.9:7;\n"
            )
            self.write_manifest(
                root,
                [
                    ("OG0000001", binary_name),
                    ("OG0000002", polytomy_name),
                ],
            )

            result = self.run_preparation(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root / binary_name).read_text(), binary)

            resolved = parse_newick((root / polytomy_name).read_text())
            self.assertEqual(
                sorted(node["label"] for node in walk(resolved) if not node["children"]),
                ["A", "B", "C", "D", "E"],
            )
            self.assertTrue(
                all(
                    len(node["children"]) == 2
                    for node in walk(resolved)
                    if node["children"]
                )
            )

            with (root / "speciesrax_gene_tree_validation.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["status"] for row in rows], ["binary", "resolved"])
            self.assertEqual(rows[1]["multifurcations_resolved"], "2")
            self.assertEqual(rows[1]["zero_length_edges_inserted"], "2")

    def test_rejects_unary_internal_nodes(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree_name = "OG0000001_ft.newick"
            (root / tree_name).write_text("((A:1):2,B:3);\n")
            self.write_manifest(root, [("OG0000001", tree_name)])

            result = self.run_preparation(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unary internal nodes", result.stderr)


if __name__ == "__main__":
    unittest.main()
