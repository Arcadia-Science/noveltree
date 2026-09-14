"""Opt-in WITCH container integration checks."""
import json

from container_helpers import CONTAINER_TESTS, INTEGRATION, invoke

pytestmark = CONTAINER_TESTS


def test_witch_nonroot_empty_store_qc_and_reuse(tmp_path):
    results = tmp_path / "results"
    (results / "alignments" / "original").mkdir(parents=True)
    args = ["--fixtures", str(INTEGRATION / "fixtures")]
    first = invoke(tmp_path, "witch.nf", results, args)
    assert "Submitted process > WITCH" in first
    alignment = results / "alignments" / "original" / "OGwitch_witch.fa"
    assert alignment.read_text().count(">") == 8
    mapping = alignment.parent / "species_protein_maps" / "OGwitch_map.link"
    assert len(mapping.read_text().splitlines()) == 8
    second = invoke(tmp_path, "witch.nf", results, args)
    assert "Stored process > WITCH" in second
    assert "Submitted process > WITCH" not in second
    # A genuinely QC-rejected task must emit a required record, not an empty MSA.
    rejected = tmp_path / "rejected"
    invoke(tmp_path, "witch.nf", rejected, args + ["--min_num_seq_per_og", "20"])
    outcome = rejected / "alignments" / "original" / "OGwitch_witch.qc_failed.json"
    assert json.loads(outcome.read_text())["reason"] == "alignment_qc"
    assert not outcome.with_name("OGwitch_witch.fa").exists()


def test_witch_with_query_sequences(tmp_path):
    fixtures = tmp_path / "fixtures"
    fixtures.mkdir()
    original = (INTEGRATION / "fixtures" / "OGwitch.fa").read_text()
    # Force one short sequence out of the median-length backbone and through
    # WITCH's HMM/query merge path, in addition to the backbone-only test.
    original = original.split(">H_p")[0] + ">H_p\nACDEFGHIKLMNPQHSTVWYACDEFGHIK\n"
    (fixtures / "OGwitch.fa").write_text(original)
    (fixtures / "OGwitch.map.link").write_text((INTEGRATION / "fixtures" / "OGwitch.map.link").read_text())
    results = tmp_path / "results"
    invoke(tmp_path, "witch.nf", results, ["--fixtures", str(fixtures)])
    assert (results / "alignments" / "original" / "OGwitch_witch.fa").read_text().count(">") == 8
