import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
class PhyloProfilesTests(unittest.TestCase):
    def test_profiles_use_coverage_without_orthofinder_matrix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            coverage = root / "coverage.txt"
            coverage.write_text(
                "SPECIES: FAMILY_COVERAGE\nSpecies-a: 1\nSpecies-b: 0.5\n"
            )
            events = root / "events.csv"
            events.write_text(
                "species_label, speciations, duplications, losses, transfers, presence, origination\n"
                "Node_Species-a_Species-b_0, 1, 2, 3, 0, 1, 0\n"
                "Species-b, 4, 5, 6, 0, 1, 0\n"
                "Species-a, 7, 8, 9, 0, 1, 0\n"
            )

            subprocess.run(
                [
                    "Rscript",
                    str(REPO / "bin" / "phylo_profiles.R"),
                    str(events),
                    str(coverage),
                    "OG1",
                ],
                cwd=root,
                check=True,
                capture_output=True,
                text=True,
            )

            with (root / "OG1_duplication_count.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                rows,
                [{"gene_family": "OG1", "Species-a": "8", "Species-b": "5"}],
            )


if __name__ == "__main__":
    unittest.main()
