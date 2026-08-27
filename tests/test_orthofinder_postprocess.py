import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


class OrthofinderPostprocessTests(unittest.TestCase):
    def test_fasta_metadata_manifest(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            species_tree = root / "species_tree"
            gene_tree = root / "gene_tree"
            species_tree.mkdir()
            gene_tree.mkdir()
            (species_tree / "OG1.fa").write_text(
                ">species-one_p1\nAAAA\nAA\n>species-two_p2\nCCC\n"
            )
            (gene_tree / "OG2.fa").write_text(
                ">species-one_p3\nA\n>species-one_p4\nTTTT\n"
            )
            output = root / "metadata.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "bin" / "summarize_og_fastas.py"),
                    "--species-tree-dir",
                    str(species_tree),
                    "--gene-tree-dir",
                    str(gene_tree),
                    "--output",
                    str(output),
                ],
                check=True,
            )

            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                rows,
                [
                    {
                        "orthogroup": "OG1",
                        "family_set": "species_tree",
                        "n_seq": "2",
                        "max_len": "6",
                        "n_species": "2",
                    },
                    {
                        "orthogroup": "OG2",
                        "family_set": "gene_tree",
                        "n_seq": "2",
                        "max_len": "4",
                        "n_species": "1",
                    },
                ],
            )

    def test_orthofinder_membership_is_not_mutated(self):
        module = (REPO / "modules/local/orthofinder_mcl.nf").read_text()
        self.assertEqual(module.count("og_tax_summary.py"), 1)
        self.assertNotIn("chimera", module.lower())
        self.assertNotIn("checkpoint", module.lower())
        self.assertNotIn("Orthogroups.post_", module)


if __name__ == "__main__":
    unittest.main()
