from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_standard_resources_match_released_process_policies() -> None:
    base = (ROOT / "conf" / "base.config").read_text()
    clipkit = (ROOT / "modules" / "local" / "clipkit.nf").read_text()
    fasttree = (ROOT / "modules" / "local" / "fasttree.nf").read_text()

    assert "check_max( 192" in base
    assert "check_max( 192.GB * task.attempt" in base
    assert "check_max( 48.h   * task.attempt" in base
    assert "check_max( 1" in base
    assert "check_max( 2.GB * task.attempt" in base

    assert "time { 6.h * task.attempt }" in clipkit
    assert "Math.min(estimated_gb, 64L)" in clipkit
    assert "fasta.size()" not in clipkit

    assert "Math.min( 16 * task.attempt" in fasttree
    assert "3.h * Math.pow(3, task.attempt - 1)" in fasttree
    assert "fasttree_large_family" not in fasttree


def test_large_run_resources_are_isolated_to_opt_in_profile() -> None:
    arcadia = (ROOT / "conf" / "arcadia.config").read_text()
    large = (ROOT / "conf" / "arcadia_large.config").read_text()
    config = (ROOT / "nextflow.config").read_text()

    assert "includeConfig 'conf/arcadia_large.config'" in config
    assert "960.GB" not in arcadia
    assert "1440.GB" not in arcadia
    assert "fasttree_large_family" not in arcadia

    assert "withLabel:process_bundle" in large
    assert "withLabel:process_orthofinder" in large
    assert "withName:SPECIESRAX" in large
    assert "withName:'FASTTREE_TIER2|FASTTREE_FALLBACK'" in large
    assert "withName:CLIPKIT" in large
    assert "memory     = 960.GB" in large
    assert "memory     = 1440.GB" in large
    assert "fasttree_large_family_min_sequences = 15000" in large


def test_one_off_fasttree_resume_auditor_is_absent() -> None:
    assert not (ROOT / "bin" / "audit_fasttree_resume.py").exists()
