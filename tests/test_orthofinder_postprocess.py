import csv
import gzip
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_script(name):
    path = REPO / "bin" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHIMERAS = load_script("flag_cross_og_chimeras.py")


class OrthofinderPostprocessTests(unittest.TestCase):
    def test_orthogroup_field_larger_than_default_csv_limit(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "Orthogroups.tsv"
            proteins = [f"speciesA_protein_{index:06d}" for index in range(8_000)]
            self.assertGreater(len(", ".join(proteins)), 131_072)
            path.write_text(
                "Orthogroup\tspeciesA\tspeciesB\n"
                f"OG0000000\t{', '.join(proteins)}\t\n"
                "OG0000001\t\tspeciesB_protein_1\n"
            )

            CHIMERAS.maximize_csv_field_size()
            protein_to_og, retained = CHIMERAS.build_protein_og_map(
                path, {"OG0000000"}
            )

            self.assertEqual(len(protein_to_og), 8_001)
            self.assertEqual(len(retained), 8_000)
            self.assertEqual(protein_to_og[proteins[-1]], "OG0000000")
            self.assertNotIn("speciesB_protein_1", retained)

    def test_only_retained_queries_are_scored_and_removed(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            blast = root / "Blast0_0.txt"
            blast.write_text(
                "a1\ta2\t0\t0\t0\t0\t0\t0\t0\t0\t1e-50\t200\n"
                "a1\tx1\t0\t0\t0\t0\t0\t0\t0\t0\t1e-40\t150\n"
                "a2\ta1\t0\t0\t0\t0\t0\t0\t0\t0\t1e-50\t200\n"
                "x1\ty1\t0\t0\t0\t0\t0\t0\t0\t0\t1e-50\t200\n"
                "x1\ta1\t0\t0\t0\t0\t0\t0\t0\t0\t1e-40\t180\n"
            )
            protein_to_og = {
                "a1": "OG_KEEP",
                "a2": "OG_KEEP",
                "x1": "OG_OTHER",
                "y1": "OG_OTHER",
            }
            retained_species = {"a1": "speciesA", "a2": "speciesA"}

            chunks = list(
                CHIMERAS.iter_blast_score_chunks(
                    root, protein_to_og, retained_species, 100, 1e-10
                )
            )
            self.assertEqual(len(chunks), 1)
            scores = chunks[0][1]
            chimeras = CHIMERAS.identify_chimeras(scores, protein_to_og, 0.5)

            self.assertEqual([record[0] for record in chimeras], ["a1"])
            self.assertNotIn("x1", scores)

            fasta_directory = root / "fastas"
            fasta_directory.mkdir()
            keep_fasta = fasta_directory / "OG_KEEP.fa"
            keep_fasta.write_text(">a1\nAAAA\n>a2\nCCCC\n")
            other_fasta = fasta_directory / "OG_OTHER.fa"
            other_fasta.write_text(">x1\nGGGG\n>y1\nTTTT\n")
            removed = CHIMERAS.remove_from_og_fastas(chimeras, fasta_directory)

            self.assertEqual(removed, 1)
            self.assertNotIn(">a1\n", keep_fasta.read_text())
            self.assertEqual(other_fasta.read_text(), ">x1\nGGGG\n>y1\nTTTT\n")

    def test_gzip_compressed_blast_results_are_scored(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            with gzip.open(root / "Blast12_7.txt.gz", "wt") as handle:
                handle.write(
                    "a1\ta2\t0\t0\t0\t0\t0\t0\t0\t0\t1e-50\t200\n"
                    "a1\tx1\t0\t0\t0\t0\t0\t0\t0\t0\t1e-40\t150\n"
                )

            chunks = list(
                CHIMERAS.iter_blast_score_chunks(
                    root,
                    {"a1": "OG_KEEP", "a2": "OG_KEEP", "x1": "OG_OTHER"},
                    {"a1": "speciesA"},
                    100,
                    1e-10,
                )
            )

            self.assertEqual([chunk[0] for chunk in chunks], [12])
            self.assertEqual(chunks[0][1]["a1"]["OG_KEEP"], 200)
            self.assertEqual(chunks[0][1]["a1"]["OG_OTHER"], 150)

    def test_fasta_metadata_manifest(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            species_tree = root / "species_tree"
            gene_tree = root / "gene_tree"
            species_tree.mkdir()
            gene_tree.mkdir()
            (species_tree / "OG1.fa").write_text(
                ">species_one_p1\nAAAA\nAA\n>species_two_p2\nCCC\n"
            )
            (gene_tree / "OG2.fa").write_text(
                ">species_one_p3\nA\n>species_one_p4\nTTTT\n"
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

    def test_post_chimera_counts_drive_final_streaming_filter(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            original_counts = root / "Orthogroups.GeneCount.tsv"
            original_counts.write_text(
                "Orthogroup\tspeciesA\tspeciesB\tTotal\n"
                "OG_KEEP\t2\t1\t3\n"
                "OG_OTHER\t1\t1\t2\n"
            )
            updated_counts = root / "post_chimera.tsv"
            chimeras = [("a1", "OG_KEEP", "OG_OTHER", "0.750")]
            CHIMERAS.write_updated_gene_counts(
                original_counts,
                updated_counts,
                chimeras,
                {"a1": "speciesA"},
            )
            with updated_counts.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["speciesA"], "1")
            self.assertEqual(rows[0]["Total"], "2")

            original_membership = root / "Orthogroups.tsv"
            original_membership.write_text(
                "Orthogroup\tspeciesA\tspeciesB\n"
                "OG_KEEP\ta1, a2\tb1\n"
                "OG_OTHER\tx1\ty1\n"
            )
            updated_membership = root / "Orthogroups.post_chimera.tsv"
            CHIMERAS.write_updated_orthogroups(
                original_membership, updated_membership, chimeras
            )
            with updated_membership.open() as handle:
                membership_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(membership_rows[0]["speciesA"], "a2")
            self.assertEqual(membership_rows[0]["speciesB"], "b1")
            self.assertEqual(membership_rows[1]["speciesA"], "x1")

            samplesheet = root / "samplesheet.csv"
            samplesheet.write_text("sample\nspeciesA\nspeciesB\n")
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "bin" / "og_tax_summary.py"),
                    str(updated_counts),
                    str(samplesheet),
                    "3",
                    "2",
                    "0.5",
                    "1",
                ],
                cwd=root,
                check=True,
                capture_output=True,
                text=True,
            )
            with (root / "spptree_core_ogs_counts.csv").open() as handle:
                species_tree_rows = list(csv.DictReader(handle))
            with (root / "genetree_core_ogs_counts.csv").open() as handle:
                gene_tree_rows = list(csv.DictReader(handle))

            self.assertEqual(species_tree_rows, [])
            self.assertEqual(gene_tree_rows, [])


if __name__ == "__main__":
    unittest.main()
