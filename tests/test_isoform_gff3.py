import csv
import gzip
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from filter_isoforms_gff3 import filter_proteome, parse_attributes
from extract_ncbi_proteome import extract
from check_samplesheet import check_samplesheet
from filter_isoforms import extract_gene_id, filter_isoforms


def feature(kind, attributes, start=1):
    return f"chr1\tNCBI\t{kind}\t{start}\t100\t.\t+\t0\t{attributes}\n"


def run_filter(tmp_path, proteins, annotation, compressed=False):
    fasta, gff = tmp_path / "input.faa", tmp_path / "annotation.gff3"
    fasta.write_text(proteins)
    if compressed:
        # Magic-byte detection also works after Nextflow stageAs removes .gz.
        with gzip.open(gff, "wt") as handle:
            handle.write(annotation)
    else:
        gff.write_text(annotation)
    return filter_proteome(fasta, gff, tmp_path / "out.fa", tmp_path / "ledger.tsv",
                           tmp_path / "qc.json", "Test-species")


def test_longest_ties_multiexon_ancestry_and_identical_distinct_genes(tmp_path):
    annotation = (
        feature("CDS", "ID=cds-a;Parent=rna-a;protein_id=XP_a.1", 1)
        + feature("CDS", "ID=cds-a;Parent=rna-a;protein_id=XP_a.1", 50)
        + feature("mRNA", "ID=rna-a;Parent=gene-a")
        + feature("gene", "ID=gene-a;Dbxref=GeneID:1;Name=SHARED")
        + feature("CDS", "Parent=gene-a;protein_id=XP_b.1")
        + feature("CDS", "Parent=gene-a;protein_id=XP_short.1")
        + feature("gene", "ID=gene-b;Dbxref=GeneID:2;Name=SHARED")
        + feature("CDS", "Parent=gene-b;protein_id=XP_other.1")
    )
    result = run_filter(tmp_path, ">XP_b.1\nAAAA\n>XP_a.1 description\nAAAA*\n"
                        ">XP_short.1\nAA\n>XP_other.1\nAAAA\n", annotation, compressed=True)
    assert result["input_proteins"] == 4
    assert result["retained"] == 2
    assert result["removed_isoforms"] == 2
    assert (tmp_path / "out.fa").read_text() == ">XP_a.1 description\nAAAA*\n>XP_other.1\nAAAA\n"
    with (tmp_path / "ledger.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["representative"] == "XP_a.1"
    assert rows[3]["gene_id"] == "GeneID:2"


def test_gff_feature_gene_id_and_escaped_multi_parent(tmp_path):
    assert parse_attributes("ID=a%2Cb;Parent=x,y") == {"ID": ("a,b",), "Parent": ("x", "y")}
    annotation = (feature("gene", "ID=a%2Cb") + feature("mRNA", "ID=t1;Parent=a%2Cb")
                  + feature("mRNA", "ID=t2;Parent=a%2Cb")
                  + feature("CDS", "Parent=t1,t2;protein_id=p"))
    assert run_filter(tmp_path, ">p\nAAA\n", annotation)["retained"] == 1
    assert "gene_feature:a,b" in (tmp_path / "ledger.tsv").read_text()


@pytest.mark.parametrize("protein,annotation,reason", [
    (">p.2\nAAA\n", feature("CDS", "protein_id=p.1;Dbxref=GeneID:1"), "unmapped"),
    (">p\nAAA\n>p\nAAAA\n", feature("CDS", "protein_id=p;Dbxref=GeneID:1"), "duplicate"),
    (">p\nAAA\n", feature("CDS", "protein_id=p;Dbxref=GeneID:1")
     + feature("CDS", "protein_id=p;Dbxref=GeneID:2"), "ambiguous"),
    (">p\nAAA\n", feature("CDS", "Parent=missing;protein_id=p;Dbxref=GeneID:1"), "Missing GFF3 parent"),
    (">p\nAAA\n", feature("mRNA", "ID=t;Parent=t")
     + feature("CDS", "Parent=t;protein_id=p"), "Cycle"),
    (">p\nAAA\n", feature("gene", "ID=g;Dbxref=GeneID:1")
     + feature("CDS", "Parent=g;protein_id=p;Dbxref=GeneID:2"), "conflicting"),
    (">p\n\n", feature("CDS", "protein_id=p;Dbxref=GeneID:1"), "empty_protein"),
])
def test_strict_failures_are_audited(tmp_path, protein, annotation, reason):
    with pytest.raises(ValueError, match=reason):
        run_filter(tmp_path, protein, annotation)
    assert not (tmp_path / "out.fa").exists()
    assert json.loads((tmp_path / "qc.json").read_text())["status"] == "failed"
    assert reason in (tmp_path / "ledger.tsv").read_text()


def test_conflicting_feature_definitions_fail(tmp_path):
    with pytest.raises(ValueError, match="Conflicting definitions"):
        run_filter(tmp_path, ">p\nAA\n", feature("gene", "ID=g;Dbxref=GeneID:1")
                   + feature("gene", "ID=g;Dbxref=GeneID:2"))


def test_ncbi_pair_selection_and_missing_annotation(tmp_path):
    root = tmp_path / "ncbi_dataset" / "data"
    assembly = root / "GCF_123.2"
    assembly.mkdir(parents=True)
    (root / "assembly_data_report.jsonl").write_text('{"accession":"GCF_123.2"}\n')
    (assembly / "protein.faa").write_text(">p\nAA\n")
    paths = (tmp_path / "p.fa", tmp_path / "g.gff", tmp_path / "meta.json")
    with pytest.raises(ValueError, match="Missing/empty"):
        extract(tmp_path, "GCF_123", *paths)
    (assembly / "genomic.gff").write_text("##gff-version 3\n")
    extract(tmp_path, "GCF_123", *paths)
    assert json.loads(paths[-1].read_text())["resolved_accession"] == "GCF_123.2"
    with pytest.raises(ValueError, match="exactly one"):
        extract(tmp_path, "GCF_123.1", *paths)


def test_samplesheet_quoted_annotations_and_old_columns(tmp_path):
    source, output = tmp_path / "samples.csv", tmp_path / "checked.csv"
    source.write_text('species,input_data,input_type,annotation_gff3,filter_isoforms\n'
                      'Test_species,input.faa,proteins,"s3://bucket/a,b.gff3.gz",yes\n')
    check_samplesheet(source, output)
    with output.open() as handle:
        row = next(csv.DictReader(handle))
    assert row["species"] == "Test-species"
    assert row["annotation_gff3"] == "s3://bucket/a,b.gff3.gz"
    source.write_text("species,input_data,input_type\nTest-species,input.fa,proteins\n")
    check_samplesheet(source, output)
    with output.open() as handle:
        assert next(csv.DictReader(handle))["annotation_gff3"] == ""


def test_annotation_with_transcriptome_is_rejected(tmp_path):
    source = tmp_path / "samples.csv"
    source.write_text("species,input_data,input_type,transdecoder,annotation_gff3\n"
                      "Test-species,input.fna,transcriptome,yes,genes.gff3\n")
    with pytest.raises(SystemExit):
        check_samplesheet(source, tmp_path / "checked.csv")


def test_legacy_brackets_gzip_and_missing_gene_warning(tmp_path, capsys):
    assert extract_gene_id("XP_1 [gene=LOC1]") == "LOC1"
    source = tmp_path / "proteins.data"
    with gzip.open(source, "wt") as handle:
        handle.write(">p1 [gene=LOC1]\nAAAA\n>p2 gene=LOC1\nAA\n>unknown\nAAA\n")
    output = tmp_path / "out.fa"
    filter_isoforms(source, output)
    assert output.read_text().count(">") == 2
    assert "1 proteins have no recognized" in capsys.readouterr().err
