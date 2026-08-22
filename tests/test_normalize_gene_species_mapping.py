import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "normalize_gene_species_mapping.py"


class NormalizeGeneSpeciesMappingTests(unittest.TestCase):
    def run_normalizer(self, root, mapping_text, tree_text=None):
        mapping = root / "family.map"
        tree = root / "species.newick"
        mapping.write_text(mapping_text)
        tree.write_text(
            tree_text
            or "((Saccharum-spontaneum:1,Oryza-sativa:1)'clade one':1,"
            "'Zea-mays':1)root;\n"
        )
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--mapping",
                str(mapping),
                "--species-tree",
                str(tree),
                "--orthogroup",
                "OG0000001",
            ],
            capture_output=True,
            text=True,
        )
        return result, mapping

    def test_repairs_malformed_species_from_gene_prefix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            gene = "Saccharum-spontaneum_Sspon.06G0004800-2B-mRNA-1"
            result, mapping = self.run_normalizer(
                root,
                f"{gene}\t{gene}\nOryza-sativa_gene_1\tOryza-sativa\n",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("2 rows validated; 1 corrected", result.stdout)
            self.assertEqual(
                mapping.read_text(),
                f"{gene}\tSaccharum-spontaneum\n"
                "Oryza-sativa_gene_1\tOryza-sativa\n",
            )

    def test_valid_mapping_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            original = "Zea-mays_gene_with_underscores\tZea-mays\n"
            result, mapping = self.run_normalizer(root, original)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("1 rows validated; 0 corrected", result.stdout)
            self.assertEqual(mapping.read_text(), original)

    def test_rejects_species_absent_from_tree(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result, _mapping = self.run_normalizer(
                root, "Unknown-species_gene1\tUnknown-species\n"
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("is not a leaf in the GeneRax species tree", result.stderr)

    def test_rejects_duplicate_genes(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result, _mapping = self.run_normalizer(
                root,
                "Oryza-sativa_gene1\tOryza-sativa\n"
                "Oryza-sativa_gene1\tOryza-sativa\n",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Duplicate gene", result.stderr)

    def test_reads_quoted_leaf_labels_and_ignores_internal_labels(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result, _mapping = self.run_normalizer(
                root,
                "Species-one_gene1\tbad\nSpecies-two_gene2\tSpecies-two\n",
                "('Species-one'[note]:1,Species-two:1)Internal-node:1;\n",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("2 rows validated; 1 corrected", result.stdout)


if __name__ == "__main__":
    unittest.main()
