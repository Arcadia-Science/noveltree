import csv
import gzip
import importlib.util
import json
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


NORMALIZE = load_script("normalize_proteome_fasta.py")
FAMILY_MAPS = load_script("build_family_maps.py")
SUBSET = load_script("subset_gene_species_map.py")


class CanonicalMappingTests(unittest.TestCase):
    def test_normalizes_gzip_ids_and_rare_amino_acids(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "input.fa.gz"
            with gzip.open(source, "wt") as handle:
                handle.write(">protein_one description\nAUOZ\n>Genus-species_already-ok\nOO\n")
            output = root / "Genus-species.fa"
            mapping = root / "Genus-species.protein_map.tsv"
            qc = root / "Genus-species.normalization_qc.json"

            NORMALIZE.normalize(
                source, "Genus-species", output, mapping, qc
            )

            self.assertEqual(
                output.read_text(),
                ">Genus-species_protein-one\nAXXZ\n"
                ">Genus-species_already-ok\nXX\n",
            )
            with mapping.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["original_protein_id"], "protein_one")
            self.assertEqual(
                rows[0]["canonical_protein_id"], "Genus-species_protein-one"
            )
            metrics = json.loads(qc.read_text())
            self.assertEqual(metrics["selenocysteine_u_replaced_with_x"], 1)
            self.assertEqual(metrics["pyrrolysine_o_replaced_with_x"], 3)

    def test_supported_uniprot_headers_preserve_clean_accessions(self):
        self.assertEqual(
            NORMALIZE.canonical_protein_id("Genus-species", "sp|P12345|ENTRY"),
            "Genus-species_P12345",
        )
        self.assertEqual(
            NORMALIZE.canonical_protein_id(
                "Genus-species", "Genus-species:Q9XYZ1"
            ),
            "Genus-species_Q9XYZ1",
        )
        self.assertEqual(
            NORMALIZE.canonical_protein_id("Genus-species", "gene:unsafe/id"),
            "Genus-species_gene-unsafe-id",
        )

    def test_rejects_canonical_identifier_collisions(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "input.fa"
            source.write_text(">a_b\nAA\n>a-b\nCC\n")
            with self.assertRaisesRegex(ValueError, "collision"):
                NORMALIZE.normalize(
                    source,
                    "Species-one",
                    root / "out.fa",
                    root / "map.tsv",
                    root / "qc.json",
                )

    def test_family_maps_and_trimmed_subsets_preserve_species(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            canonical = root / "Species-one.protein_map.tsv"
            canonical.write_text(
                "species\toriginal_protein_id\tcanonical_protein_id\n"
                "Species-one\ta\tSpecies-one_a\n"
                "Species-one\tb\tSpecies-one_b\n"
            )
            species_tree = root / "species_tree"
            gene_tree = root / "gene_tree"
            family_maps = root / "family_maps"
            species_tree.mkdir()
            gene_tree.mkdir()
            (species_tree / "OG1.fa").write_text(
                ">Species-one_a\nAA\n>Species-one_b\nCC\n"
            )
            previous = Path.cwd()
            try:
                import os

                os.chdir(root)
                lookup = FAMILY_MAPS.read_canonical_maps("*.protein_map.tsv")
                FAMILY_MAPS.write_family_maps(
                    lookup, [species_tree, gene_tree], family_maps
                )
            finally:
                os.chdir(previous)

            family_map = family_maps / "OG1_map.link"
            self.assertEqual(
                family_map.read_text(),
                "Species-one_a\tSpecies-one\nSpecies-one_b\tSpecies-one\n",
            )
            trimmed = root / "trimmed.fa"
            trimmed.write_text(">Species-one_b\nCC\n")
            subset = root / "trimmed.map.link"
            SUBSET.subset(family_map, trimmed, subset)
            self.assertEqual(subset.read_text(), "Species-one_b\tSpecies-one\n")


if __name__ == "__main__":
    unittest.main()
