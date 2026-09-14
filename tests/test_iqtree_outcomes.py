import json
from pathlib import Path
import shutil
import subprocess

import pytest

RUNNER = Path(__file__).resolve().parents[1] / "bin" / "run_iqtree_with_fallback.sh"
pytestmark = pytest.mark.skipif(not shutil.which("timeout"), reason="GNU timeout required")


def run(tmp_path, command, budget="5"):
    return subprocess.run(["bash", str(RUNNER), "OG1", budget, "1", "source.treefile",
                           "final.newick", "fallback.json", "--", "bash", "-c", command],
                          cwd=tmp_path, capture_output=True, text=True)


def test_success(tmp_path):
    assert run(tmp_path, "printf '(A,B);' > source.treefile").returncode == 0
    assert (tmp_path / "final.newick").read_text() == "(A,B);"
    assert not (tmp_path / "fallback.json").exists()


@pytest.mark.parametrize("command,reason,exit_status", [
    ("exit 1", "tool_exit", 1),
    ("exit 139", "tool_exit", 139),
    ("true", "missing_tree", 0),
    ("printf '(partial' > source.treefile; sleep 10", "timeout", 124),
])
def test_handled_failure_never_publishes_partial_tree(tmp_path, command, reason, exit_status):
    assert run(tmp_path, command, budget="1").returncode == 0
    record = json.loads((tmp_path / "fallback.json").read_text())
    assert record["reason"] == reason
    assert record["exit_status"] == exit_status
    assert not (tmp_path / "final.newick").exists()


def test_missing_executable_remains_real_failure(tmp_path):
    result = subprocess.run(["bash", str(RUNNER), "OG1", "1", "1", "source",
                             "final", "fallback.json", "--", "nonexistent-iqtree-command"],
                            cwd=tmp_path, capture_output=True)
    assert result.returncode == 127
    assert not (tmp_path / "fallback.json").exists()
