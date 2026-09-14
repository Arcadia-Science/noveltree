"""Opt-in real Nextflow dataflow tests; only expensive scientific executables are doubled."""
import json
import os
from pathlib import Path
import subprocess

import pytest

REPO = Path(__file__).resolve().parents[1]
INTEGRATION = REPO / "tests" / "integration"
pytestmark = pytest.mark.skipif(os.environ.get("RUN_NEXTFLOW_TESTS") != "1",
                               reason="set RUN_NEXTFLOW_TESTS=1 (requires Nextflow + Java + GNU timeout)")

FASTA = ">A_p\nACDEFGHIKLMNPQRSTVWY\n>B_p\nACDEFGHIKLMNPQRSTVWY\n>C_p\nACDEFGHIKLMNPQRSTVWY\n>D_p\nACDEFGHIKLMNPQRSTVWY\n"
TREE = "((A_p:0.1,B_p:0.1):0.1,(C_p:0.1,D_p:0.1):0.1);\n"


def test_streaming_fallback_stored_outputs_and_external_failure(tmp_path):
    inputs, results, events = [tmp_path / name for name in ("inputs", "results", "events")]
    for directory in (inputs, results, events):
        directory.mkdir()
    alignments = results / "alignments" / "original"
    maps = alignments / "species_protein_maps"
    maps.mkdir(parents=True)
    original = results / "gene_family_trees" / "original"
    original.mkdir(parents=True)
    ogs = ["OGfail", "OGslow", "OGtimeout", "OGexternal", "OGstored", "OGcached", "OGlarge"]
    for og in ogs:
        (inputs / f"{og}.fa").write_text(FASTA)
        mapping = "".join(f"{species}_p\t{species}\n" for species in "ABCD")
        (inputs / f"{og}.map.link").write_text(mapping)
        (alignments / f"{og}_famsa.fa").write_text(FASTA)
        (maps / f"{og}_map.link").write_text(mapping)
    # Old successful stored files, deliberately WITHOUT the new JSON sidecars.
    (original / "OGstored_famsa_ft.newick").write_text(TREE)
    (original / "OGcached_famsa_iqt.newick").write_text(TREE)
    command = ["nextflow", "run", str(INTEGRATION / "trees.nf"), "-ansi-log", "false",
               "--fixtures", str(inputs), "--outdir", str(results), "-work-dir", str(tmp_path / "work")]
    result = subprocess.run(command, cwd=tmp_path, env=dict(os.environ, FAKE_EVENT_DIR=str(events)),
                            capture_output=True, text=True, timeout=120)
    (tmp_path / "run.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert float((events / "OGfail.ft_start").read_text()) < float((events / "OGslow.iq_done").read_text())
    for og in ("OGfail", "OGtimeout", "OGexternal", "OGlarge"):
        assert len((events / f"{og}.ft_start").read_text().splitlines()) == 1
    for og in ("OGstored", "OGcached"):
        assert not (events / f"{og}.iq_start").exists()
        assert not (events / f"{og}.ft_start").exists()
    for og in ogs:
        assert len(list(original.glob(f"{og}_*.newick"))) == 1
        assert (results / "gene_family_trees" / "reconciliation_ready" / f"{og}_reconciliation_ready.newick").exists()
    assert json.loads((original / "OGtimeout_famsa_iqt.fallback.json").read_text())["reason"] == "timeout"
    assert not (original / "OGtimeout_famsa_iqt.newick").exists()
