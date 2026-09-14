"""Helpers for optional container integration tests."""
import os
from pathlib import Path
import subprocess

import pytest

INTEGRATION = Path(__file__).resolve().parent / "integration"
CONTAINER_TESTS = pytest.mark.skipif(
    os.environ.get("RUN_DOCKER_TESTS") != "1",
    reason="set RUN_DOCKER_TESTS=1 (Nextflow, Docker and module images required)",
)


def invoke(tmp_path, workflow, results, arguments=()):
    result = subprocess.run(
        ["nextflow", "run", str(INTEGRATION / workflow), "-c", str(INTEGRATION / "docker.config"),
         "-ansi-log", "false", "--outdir", str(results), "-work-dir", str(tmp_path / "work"),
         *arguments], cwd=tmp_path, capture_output=True, text=True, timeout=240)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr
