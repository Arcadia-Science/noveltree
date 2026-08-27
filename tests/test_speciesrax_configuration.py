from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_undated_dl_speciesrax_does_not_enable_species_tree_pruning():
    config = (REPO_ROOT / "conf" / "modules.config").read_text()
    speciesrax_blocks = re.findall(
        r"withName:\s*'SPECIESRAX'\s*\{(.*?)\n\s*\}",
        config,
        flags=re.DOTALL,
    )

    assert len(speciesrax_blocks) == 2
    for block in speciesrax_blocks:
        assert "--rec-model UndatedDL" in block
        assert "--prune-species-tree" not in block


def test_speciesrax_rejects_the_known_broken_generax_combination():
    module = (REPO_ROOT / "modules" / "local" / "speciesrax.nf").read_text()

    assert "args.contains('--prune-species-tree')" in module
    assert "args.contains('UndatedDL')" in module


def test_speciesrax_preserves_an_explicit_reference_root():
    config = (REPO_ROOT / "conf" / "modules.config").read_text()
    module = (REPO_ROOT / "modules" / "local" / "speciesrax.nf").read_text()
    main = (REPO_ROOT / "main.nf").read_text()

    assert "REROOT" not in config
    assert config.count("--si-strategy SKIP") == 2
    assert "--species-tree MiniNJ" in module
    assert "--do-not-reconcile" in module
    assert "root_species_tree_from_reference.py" in module
    assert "--species-tree rooted_mininj_species_tree.newick" in module
    assert "BUILD_REFERENCE_CHRONOGRAM" in main
    assert "Channel.value([])" in main
    assert "no_reference_chronogram.sentinel" not in main
    assert "path reference_chronogram" in module
    assert "Channel.value([])" in main


def test_reference_root_transfer_validates_the_chronogram():
    script = (REPO_ROOT / "bin" / "root_species_tree_from_reference.py").read_text()
    calibration = (
        REPO_ROOT / "bin" / "zoogle" / "time_calibrate_species_tree.R"
    ).read_text()

    assert "must have exactly two children at its encoded root" in script
    assert "is not ultrametric at its encoded root" in script
    assert "root_split_validated" in script
    assert "SpeciesRax root split does not match the reference chronogram" in calibration
