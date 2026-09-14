"""Opt-in proteome preprocessing container integration checks."""
import csv
import gzip
import json

from container_helpers import CONTAINER_TESTS, invoke

pytestmark = CONTAINER_TESTS


def test_gff_preprocessing_keeps_identical_distinct_genes(tmp_path):
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
    proteins, annotation = tmp_path / "proteins.faa.gz", tmp_path / "genes,annotation.gff3.gz"
    with gzip.open(proteins, "wt") as handle:
        handle.write(f">p1\n{sequence[:30]}\n>p2\n{sequence}\n>p3\n{sequence}\n")
    with gzip.open(annotation, "wt") as handle:
        for protein, gene in [("p1", 1), ("p2", 1), ("p3", 2)]:
            handle.write(f"chr1\tNCBI\tCDS\t1\t120\t.\t+\t0\tprotein_id={protein};Dbxref=GeneID:{gene}\n")
    legacy = tmp_path / "legacy.fa"
    legacy.write_text(f">a [gene=g1]\n{sequence}\n>b [gene=g2]\n{sequence}\n")
    samplesheet = tmp_path / "samples.csv"
    with samplesheet.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", "input_data", "input_type", "filter_isoforms", "annotation_gff3"])
        writer.writerow(["Annotated-species", str(proteins), "proteins", "yes", str(annotation)])
        writer.writerow(["Legacy-species", str(legacy), "proteins", "yes", ""])
    results = tmp_path / "results"
    invoke(tmp_path, "inputs.nf", results, ["--input", str(samplesheet)])
    final = results / "preprocessing" / "final_proteomes"
    assert (final / "Annotated-species.fa").read_text().count(">") == 2
    assert (final / "Legacy-species.fa").read_text().count(">") == 1
    qc = json.loads((results / "preprocessing" / "isofiltered" / "Annotated-species_isoform_qc.json").read_text())
    assert qc["retained"] == 2
    assert qc["removed_isoforms"] == 1
