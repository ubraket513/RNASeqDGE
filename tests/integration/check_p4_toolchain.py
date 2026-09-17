"""Offline staging/preflight safety checks; no third-party Python packages."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import subprocess

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("toolchain", ROOT / "tools/toolchain.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ToolchainTests(unittest.TestCase):
    def test_installed_records_require_exact_set_build_and_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            (prefix / "conda-meta").mkdir()
            record = {"name": "pkg", "version": "1", "build": "a", "sha256": "a" * 64}
            file = prefix / "conda-meta/pkg.json"
            file.write_text(json.dumps(record))
            module.check_records(prefix, [record])
            for key, value in (("version", "2"), ("build", "b"), ("sha256", "b" * 64)):
                file.write_text(json.dumps(dict(record, **{key: value})))
                with self.assertRaises(ValueError):
                    module.check_records(prefix, [record])
            file.write_text(json.dumps(record))
            (prefix / "conda-meta/extra.json").write_text(json.dumps(dict(record, name="extra")))
            with self.assertRaises(ValueError):
                module.check_records(prefix, [record])

    def test_lock_consistency_checks_url_and_checksum(self):
        record = {"url": "https://example.org/pkg.conda", "md5": "a" * 32}
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "lock.txt"
            lock.write_text("@EXPLICIT\n" + record["url"] + "#" + record["md5"] + "\n")
            module.check_lock(lock, [record])
            lock.write_text("@EXPLICIT\n" + record["url"] + "#" + "b" * 32 + "\n")
            with self.assertRaises(ValueError):
                module.check_lock(lock, [record])

    def test_existing_destination_and_dangling_link_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError):
                module.new_destination(root)
            link = root / "link"
            link.symlink_to(root / "missing")
            with self.assertRaises(ValueError):
                module.new_destination(link)
            self.assertEqual(module.new_destination(root / "new"), root / "new")

    def test_failed_staging_propagates_without_verifying_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = SimpleNamespace(prefix=Path(tmp) / "space ; literal", kind="alignment", mamba="mamba")
            with patch.object(module.shutil, "which", return_value="/tools/mamba"), \
                 patch.object(module.subprocess, "run", side_effect=subprocess.CalledProcessError(7, "mamba")) as run, \
                 patch.object(module, "check_records") as verified:
                with self.assertRaises(subprocess.CalledProcessError):
                    module.stage(args)
                argv = run.call_args.args[0]
                self.assertEqual(argv[argv.index("--prefix") + 1], str(args.prefix))
                self.assertNotIn("shell", run.call_args.kwargs)
                verified.assert_not_called()

    def test_r_environment_and_auxiliary_threads_are_constrained(self):
        with patch.dict(module.os.environ, {"R_LIBS_USER": "/untrusted", "LD_LIBRARY_PATH": "/wrong", "OMP_NUM_THREADS": "999"}):
            env = module.environment()
        self.assertEqual(env["R_LIBS_USER"], module.os.devnull)
        self.assertNotIn("LD_LIBRARY_PATH", env)
        self.assertEqual(env["OMP_NUM_THREADS"], "1")


if __name__ == "__main__":
    unittest.main()
