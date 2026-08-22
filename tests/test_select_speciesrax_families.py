import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


class SelectSpeciesRaxFamiliesTests(unittest.TestCase):
    def write_inputs(self, root, families, leaf_overrides=None):
        leaf_overrides = leaf_overrides or {}
        manifest = root / "speciesrax_inputs_00000000.tsv"
        validation = root / "speciesrax_gene_tree_validation.tsv"
        with manifest.open("w", newline="") as manifest_handle, validation.open(
            "w", newline=""
        ) as validation_handle:
            manifest_writer = csv.writer(
                manifest_handle, delimiter="\t", lineterminator="\n"
            )
            manifest_writer.writerow(["orthogroup", "gene_tree", "mapping"])
            validation_writer = csv.DictWriter(
                validation_handle,
                fieldnames=["orthogroup", "leaves"],
                delimiter="\t",
                lineterminator="\n",
            )
            validation_writer.writeheader()
            for orthogroup, species in families.items():
                tree = f"{orthogroup}.newick"
                mapping = f"{orthogroup}.map"
                (root / tree).write_text("(A:1,B:1);\n")
                with (root / mapping).open("w") as mapping_handle:
                    for index, name in enumerate(species):
                        mapping_handle.write(
                            f"{name}_{orthogroup}_gene{index} {name}\n"
                        )
                manifest_writer.writerow([orthogroup, tree, mapping])
                validation_writer.writerow(
                    {
                        "orthogroup": orthogroup,
                        "leaves": leaf_overrides.get(orthogroup, len(species)),
                    }
                )

    def run_selection(self, root):
        (root / "expected_species.txt").write_text("A\nB\nC\nD\n")
        return subprocess.run(
            [
                sys.executable,
                str(REPO / "bin" / "select_speciesrax_families.py"),
                "--manifest-glob",
                "speciesrax_inputs_*.tsv",
                "--validation-report",
                "speciesrax_gene_tree_validation.tsv",
                "--expected-species-file",
                "expected_species.txt",
                "--min-species-occupancy",
                "0.5",
                "--max-mean-copies",
                "2",
                "--max-copies-per-species",
                "3",
                "--max-total-leaves-factor",
                "2",
                "--selected-manifest",
                "speciesrax_selected_families.tsv",
                "--report",
                "speciesrax_family_selection.tsv",
                "--species-coverage-report",
                "speciesrax_selected_species_coverage.tsv",
            ],
            cwd=root,
            capture_output=True,
            text=True,
        )

    def test_applies_each_filter_and_writes_auditable_outputs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.write_inputs(
                root,
                {
                    "OG_KEEP": ["A", "A", "B", "C", "D"],
                    "OG_LOW_OCCUPANCY": ["A"],
                    "OG_HIGH_MEAN": ["A", "A", "A", "B", "B", "B"],
                    "OG_HIGH_SPECIES_COPY": ["A", "A", "A", "A", "B", "C", "D"],
                    "OG_TOO_MANY_LEAVES": [
                        "A",
                        "A",
                        "A",
                        "B",
                        "B",
                        "B",
                        "C",
                        "C",
                        "D",
                    ],
                },
            )

            result = self.run_selection(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("1/5 retained", result.stdout)

            with (root / "speciesrax_selected_families.tsv").open() as handle:
                selected = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["orthogroup"] for row in selected], ["OG_KEEP"])

            with (root / "speciesrax_family_selection.tsv").open() as handle:
                report = {
                    row["orthogroup"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(report["OG_KEEP"]["selected"], "true")
            self.assertEqual(
                report["OG_LOW_OCCUPANCY"]["exclusion_reasons"],
                "species_occupancy_below_minimum",
            )
            self.assertIn(
                "mean_copies_above_maximum",
                report["OG_HIGH_MEAN"]["exclusion_reasons"],
            )
            self.assertEqual(
                report["OG_HIGH_SPECIES_COPY"]["exclusion_reasons"],
                "species_copy_count_above_maximum",
            )
            self.assertIn(
                "total_leaves_above_maximum",
                report["OG_TOO_MANY_LEAVES"]["exclusion_reasons"],
            )

            with (root / "speciesrax_selected_species_coverage.tsv").open() as handle:
                coverage = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(coverage), 4)
            self.assertTrue(all(int(row["selected_families"]) == 1 for row in coverage))

    def test_rejects_tree_mapping_leaf_count_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.write_inputs(
                root,
                {"OG_BAD": ["A", "B", "C"]},
                leaf_overrides={"OG_BAD": 4},
            )
            result = self.run_selection(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Tree/mapping leaf-count mismatch", result.stderr)

    def test_requires_the_complete_run_level_species_set(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.write_inputs(root, {"OG_INCOMPLETE": ["A", "B", "C"]})
            result = self.run_selection(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "input species absent from all staged families: D", result.stderr
            )

    def test_repairs_malformed_species_from_gene_prefix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.write_inputs(root, {"OG_REPAIR": ["A", "B", "C", "D"]})
            mapping = root / "OG_REPAIR.map"
            lines = mapping.read_text().splitlines()
            gene = lines[0].split()[0]
            lines[0] = f"{gene}\t{gene}_cds"
            mapping.write_text("\n".join(lines) + "\n")

            result = self.run_selection(root)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(
                "1 rows corrected across 1 task-local mapping files", result.stdout
            )
            self.assertEqual(mapping.read_text().splitlines()[0], f"{gene}\tA")


if __name__ == "__main__":
    unittest.main()
