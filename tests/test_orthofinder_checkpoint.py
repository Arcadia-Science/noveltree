import argparse
import importlib.util
import io
import json
import os
import shutil
import tempfile
import unittest
from contextlib import contextmanager, redirect_stdout
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_script(name):
    path = REPO / "bin" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHECKPOINT = load_script("orthofinder_checkpoint.py")


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def make_results(root):
    results = root / "OrthoFinder" / "Results_Inflation_1.5"
    orthogroups = results / "Orthogroups"
    sequences = results / "Orthogroup_Sequences"
    orthogroups.mkdir(parents=True)
    sequences.mkdir()
    (orthogroups / "Orthogroups.tsv").write_text(
        "Orthogroup\tspeciesA\nOG0000000\tproteinA\n"
    )
    (orthogroups / "Orthogroups.GeneCount.tsv").write_text(
        "Orthogroup\tspeciesA\tTotal\nOG0000000\t1\t1\n"
    )
    (sequences / "OG0000000.fa").write_text(">proteinA\nAAAA\n")
    return results


class OrthofinderCheckpointTests(unittest.TestCase):
    def test_checkpoint_is_saved_before_custom_postprocessing(self):
        module = (
            REPO / "modules" / "local" / "orthofinder_mcl.nf"
        ).read_text()

        save_call = module.index("orthofinder_checkpoint.py save")
        chimera_call = module.index("flag_cross_og_chimeras.py")
        self.assertLess(save_call, chimera_call)
        self.assertIn("orthofinder_checkpoint.py restore", module)

    def test_fingerprint_changes_with_content_options_and_bundle_metadata(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "SpeciesIDs.txt"
            second = root / "Species0.fa"
            bundle = root / "BlastBundle_0.tar"
            first.write_text("0: speciesA.fa\n")
            second.write_text(">proteinA\nAAAA\n")
            bundle.write_bytes(b"blast-data")

            def fingerprint(inflation="1.5"):
                arguments = argparse.Namespace(
                    inputs=[str(second), str(first)],
                    metadata_only_inputs=[str(bundle)],
                    orthofinder_version="2.5.4",
                    inflation=inflation,
                    orthofinder_options="-b -I -M msa -X -os -z",
                    extra_args="",
                )
                output = io.StringIO()
                with redirect_stdout(output):
                    CHECKPOINT.fingerprint_inputs(arguments)
                return output.getvalue().strip()

            original = fingerprint()
            self.assertEqual(original, fingerprint())
            self.assertNotEqual(original, fingerprint("2.0"))
            second.write_text(">proteinA\nCCCC\n")
            self.assertNotEqual(original, fingerprint())
            second.write_text(">proteinA\nAAAA\n")
            bundle.write_bytes(b"longer-blast-data")
            self.assertNotEqual(original, fingerprint())

    def test_local_checkpoint_round_trip_and_cleanup(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            task = root / "task"
            task.mkdir()
            checkpoint_root = root / "checkpoints"
            results = make_results(task)
            receipt = task / "orthofinder_checkpoint_receipt.json"
            arguments = argparse.Namespace(
                root=str(checkpoint_root),
                fingerprint="abc123",
                results_dir=str(results.relative_to(task)),
                receipt=str(receipt),
            )

            with working_directory(task):
                self.assertEqual(CHECKPOINT.save_checkpoint(arguments), 0)
            saved_receipt = json.loads(receipt.read_text())
            self.assertTrue(Path(saved_receipt["archive"]).is_file())
            self.assertTrue(Path(saved_receipt["manifest"]).is_file())

            # A failed upload of the tiny completion manifest must not make a
            # successfully uploaded, content-addressed archive unusable.
            Path(saved_receipt["manifest"]).unlink()
            shutil.rmtree(task / "OrthoFinder")
            with working_directory(task):
                self.assertEqual(CHECKPOINT.restore_checkpoint(arguments), 0)
            self.assertEqual(
                (results / "Orthogroup_Sequences" / "OG0000000.fa").read_text(),
                ">proteinA\nAAAA\n",
            )

            cleanup = argparse.Namespace(receipt=str(receipt))
            self.assertEqual(CHECKPOINT.cleanup_checkpoint(cleanup), 0)
            self.assertFalse(Path(saved_receipt["archive"]).exists())
            self.assertFalse(Path(saved_receipt["manifest"]).exists())

    def test_missing_checkpoint_requests_fresh_clustering(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            arguments = argparse.Namespace(
                root=str(root / "missing"),
                fingerprint="abc123",
                results_dir=str(root / "OrthoFinder" / "Results_Inflation_1.5"),
                receipt=str(root / "receipt.json"),
            )
            self.assertEqual(
                CHECKPOINT.restore_checkpoint(arguments), CHECKPOINT.NOT_RESTORABLE
            )


if __name__ == "__main__":
    unittest.main()
