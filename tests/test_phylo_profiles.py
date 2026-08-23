import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
class PhyloProfilesTests(unittest.TestCase):
    def run_profiles(self, root, coverage_text, events_text):
        coverage = root / "coverage.txt"
        coverage.write_text(coverage_text)
        events = root / "events.csv"
        events.write_text(events_text)
        return subprocess.run(
            [
                "Rscript",
                str(REPO / "bin" / "phylo_profiles.R"),
                str(events),
                str(coverage),
                "OG1",
            ],
            cwd=root,
            capture_output=True,
            text=True,
        )

    def test_profiles_use_coverage_without_orthofinder_matrix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result = self.run_profiles(
                root,
                "SPECIES: FAMILY_COVERAGE\nSpecies-a: 1\nSpecies-b: 0.5\n",
                "species_label, speciations, duplications, losses, transfers, presence, origination\n"
                "Node_Species-a_Species-b_0, 1, 2, 3, 0, 1, 0\n"
                "Species-b, 4, 5, 6, 0, 1, 0\n"
                "Species-a, 7, 8, 9, 0, 1, 0\n",
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            with (root / "OG1_duplication_count.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                rows,
                [{"gene_family": "OG1", "Species-a": "8", "Species-b": "5"}],
            )

    def test_zero_fills_only_omitted_zero_coverage_species(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result = self.run_profiles(
                root,
                "SPECIES: FAMILY_COVERAGE\n"
                "Species-a: 1\nSpecies-b: 0\nSpecies-c: 0\n",
                "species_label, speciations, duplications, losses, transfers, presence, origination\n"
                "Species-a, 1, 2, 0, 0, 1, 0\n"
                "Species-b, 0, 0, 3, 0, 0, 0\n",
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            with (root / "OG1_loss_count.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                rows,
                [
                    {
                        "gene_family": "OG1",
                        "Species-a": "0",
                        "Species-b": "3",
                        "Species-c": "0",
                    }
                ],
            )

    def test_rejects_omitted_positive_coverage_species(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result = self.run_profiles(
                root,
                "SPECIES: FAMILY_COVERAGE\nSpecies-a: 1\nSpecies-b: 1\n",
                "species_label, speciations, duplications, losses, transfers, presence, origination\n"
                "Species-a, 1, 0, 0, 0, 1, 0\n",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Covered species absent from GeneRax event table: Species-b",
                result.stderr,
            )


if __name__ == "__main__":
    unittest.main()
