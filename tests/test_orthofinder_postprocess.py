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

    def test_speciesrax_constraints_route_rejected_families_to_gene_trees(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            counts = root / "Orthogroups.GeneCount.tsv"
            counts.write_text(
                "Orthogroup\tA\tB\tC\tD\tTotal\n"
                "OG_CORE\t1\t1\t1\t1\t4\n"
                "OG_OCC_BOUNDARY\t2\t2\t0\t0\t4\n"
                "OG_MAX_COPY_FAIL\t9\t1\t1\t1\t12\n"
            )
            samplesheet = root / "samplesheet.csv"
            samplesheet.write_text(
                "species,input_data,input_type\n"
                "A,a.fa,proteins\nB,b.fa,proteins\n"
                "C,c.fa,proteins\nD,d.fa,proteins\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "bin/og_tax_summary.py"),
                    str(counts),
                    str(samplesheet),
                    "4",
                    "2",
                    "0.5",
                    "10",
                    "0.5",
                    "4",
                    "8",
                    "4",
                ],
                cwd=root,
                check=True,
            )
            with (root / "spptree_core_ogs_counts.csv").open() as handle:
                species_tree = [row["orthogroup"] for row in csv.DictReader(handle)]
            with (root / "genetree_core_ogs_counts.csv").open() as handle:
                gene_tree = [row["orthogroup"] for row in csv.DictReader(handle)]
            self.assertEqual(species_tree, ["OG_CORE", "OG_OCC_BOUNDARY"])
            self.assertEqual(gene_tree, ["OG_MAX_COPY_FAIL"])

            with (root / "speciesrax_family_selection.tsv").open() as handle:
                report = {
                    row["orthogroup"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(report["OG_OCC_BOUNDARY"]["selected"], "true")
            self.assertEqual(
                report["OG_MAX_COPY_FAIL"]["exclusion_reasons"],
                "species_copy_count_above_maximum",
            )


if __name__ == "__main__":
    unittest.main()
