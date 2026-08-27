from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "audit_fasttree_resume",
    ROOT / "bin" / "audit_fasttree_resume.py",
)
assert SPEC is not None and SPEC.loader is not None
AUDIT = module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


def write_metadata(path: Path) -> None:
    path.write_text(
        "orthogroup\tfamily_set\tn_seq\tmax_len\tn_species\n"
        "OG0000000\tgene_tree\t2000\t500\t20\n"
        "OG0000001\tgene_tree\t1001\t500\t20\n"
        "OG0000002\tgene_tree\t1000\t500\t20\n"
        "OG0000003\tspecies_tree\t3000\t500\t20\n"
    )


def test_resume_audit_reports_only_nonempty_stored_outputs(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    store = tmp_path / "trees"
    store.mkdir()
    write_metadata(metadata)
    (store / "OG0000000_famsa_clipkit_ft.newick").write_text("(A,B);\n")
    (store / "OG0000001_famsa_clipkit_ft.newick").touch()

    result = AUDIT.audit_fasttree_resume(str(metadata), str(store))

    assert result.expected == (
        "OG0000000_famsa_clipkit_ft.newick",
        "OG0000001_famsa_clipkit_ft.newick",
    )
    assert result.stored == ("OG0000000_famsa_clipkit_ft.newick",)
    assert result.missing == ("OG0000001_famsa_clipkit_ft.newick",)
    assert len(result.fingerprint) == 64


def test_resume_audit_enforces_exact_expected_state(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    store = tmp_path / "trees"
    store.mkdir()
    write_metadata(metadata)
    (store / "OG0000000_famsa_clipkit_ft.newick").write_text("(A,B);\n")

    exit_code = AUDIT.main(
        [
            "--metadata",
            str(metadata),
            "--tree-store",
            str(store),
            "--expect-candidates",
            "2",
            "--expect-stored",
            "1",
            "--expect-missing",
            "OG0000001",
        ]
    )

    assert exit_code == 0


def test_snapshot_verification_ignores_new_outputs_but_detects_changes(
    tmp_path: Path,
) -> None:
    metadata = tmp_path / "metadata.tsv"
    store = tmp_path / "trees"
    snapshot = tmp_path / "snapshot.json"
    store.mkdir()
    write_metadata(metadata)
    first = store / "OG0000000_famsa_clipkit_ft.newick"
    second = store / "OG0000001_famsa_clipkit_ft.newick"
    first.write_text("(A,B);\n")

    before = AUDIT.audit_fasttree_resume(str(metadata), str(store))
    AUDIT.write_snapshot(str(snapshot), before)
    second.write_text("(C,D);\n")
    after = AUDIT.audit_fasttree_resume(str(metadata), str(store))

    assert AUDIT.verify_snapshot(str(snapshot), after) == ()

    first.write_text("(A,C);\n")
    changed = AUDIT.audit_fasttree_resume(str(metadata), str(store))
    assert AUDIT.verify_snapshot(str(snapshot), changed) == (
        "OG0000000_famsa_clipkit_ft.newick",
    )


def test_large_family_scaling_is_portable_and_not_arcadia_name_specific() -> None:
    arcadia = (ROOT / "conf" / "arcadia.config").read_text()
    module = (ROOT / "modules" / "local" / "fasttree.nf").read_text()

    assert "REMAINING_GENE_FAMILIES:FASTTREE_TIER2" not in arcadia
    assert "fasttree_large_family_min_sequences" in module
    assert "fasttree_large_family_memory_multiplier" in module
    assert "fasttree_large_family_retry_cpu_step" in module
    assert "params.awsqueue_ondemand" not in module
