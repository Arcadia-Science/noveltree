import csv
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


class SpeciesRaxBundleTests(unittest.TestCase):
    def test_bundle_contains_validated_manifest_and_inputs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for og in ("OG0000001", "OG0000002"):
                (root / f"{og}_einsi_clipkit_iqt.newick").write_text("(sp_p:1);\n")
                (root / f"{og}_einsi_map.link").write_text("sp_p sp\n")

            archive = root / "speciesrax_inputs_test.tar"
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "bin" / "bundle_speciesrax_inputs.py"),
                    "--input-dir",
                    str(root),
                    "--output",
                    str(archive),
                    "--manifest-name",
                    "speciesrax_inputs_test.tsv",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with tarfile.open(archive) as handle:
                self.assertEqual(
                    handle.getnames(),
                    [
                        "speciesrax_inputs_test.tsv",
                        "OG0000001_einsi_clipkit_iqt.newick",
                        "OG0000001_einsi_map.link",
                        "OG0000002_einsi_clipkit_iqt.newick",
                        "OG0000002_einsi_map.link",
                    ],
                )
                manifest = handle.extractfile("speciesrax_inputs_test.tsv")
                rows = list(
                    csv.DictReader(
                        (line.decode() for line in manifest), delimiter="\t"
                    )
                )
            self.assertEqual(
                [row["orthogroup"] for row in rows],
                ["OG0000001", "OG0000002"],
            )

    def test_bundle_rejects_incomplete_family(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            (root / "OG0000001_einsi_clipkit_iqt.newick").write_text("(sp_p:1);\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "bin" / "bundle_speciesrax_inputs.py"),
                    "--input-dir",
                    str(root),
                    "--output",
                    str(root / "out.tar"),
                    "--manifest-name",
                    "manifest.tsv",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("lack exactly one gene tree", result.stderr)


if __name__ == "__main__":
    unittest.main()
