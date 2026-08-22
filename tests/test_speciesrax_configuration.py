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
