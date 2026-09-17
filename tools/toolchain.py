#!/usr/bin/env python3
"""Explicit online environment staging and read-only offline toolchain preflight."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def packages(kind):
    return json.loads((ROOT / f"config/{kind}-linux-64.provenance.json").read_text())["packages"]


def check_records(prefix, expected):
    actual = [json.loads(p.read_text()) for p in (prefix / "conda-meta").glob("*.json")]
    fields = ("name", "version", "build", "sha256")
    def records(rows):
        result = {}
        for row in rows:
            if not all(row.get(key) for key in fields) or row["name"] in result:
                raise ValueError(f"{prefix}: incomplete or duplicate package record")
            result[row["name"]] = tuple(row[key] for key in fields)
        return result
    wanted, found = records(expected), records(actual)
    if wanted != found:
        changed = sorted(k for k in wanted.keys() | found.keys() if wanted.get(k) != found.get(k))
        raise ValueError(f"{prefix}: package lock mismatch: {', '.join(changed)}")


def check_lock(path, expected):
    lines = [line.strip() for line in path.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    wanted = [row["url"] + "#" + row["md5"] for row in expected]
    if not lines or lines[0] != "@EXPLICIT" or sorted(lines[1:]) != sorted(wanted):
        raise ValueError(f"{path}: explicit URLs/checksums do not match provenance")


def new_destination(path):
    path = Path(os.path.abspath(path))
    if path.exists() or path.is_symlink():
        raise ValueError(f"destination already exists: {path}")
    if not path.parent.is_dir():
        raise ValueError(f"destination parent must exist: {path.parent}")
    return path


def environment():
    env = dict(os.environ)
    # Keep hidden library paths and user R profiles out of the pinned preflight.
    for key in ("LD_LIBRARY_PATH", "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"):
        env.pop(key, None)
    # An unset R_LIBS_USER enables R's default ~/R library; /dev/null is never
    # a library directory, while the selected R installation keeps .Library.
    env["R_LIBS_USER"] = os.devnull
    env["R_LIBS_SITE"] = os.devnull
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"
    return env


def stage(args):
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "amd64"):
        raise ValueError("locks support Linux x86_64 only")
    prefix = new_destination(args.prefix)
    lock = ROOT / f"config/{args.kind}-linux-64.explicit.txt"
    check_lock(lock, packages(args.kind))
    manager = shutil.which(args.mamba)
    if manager is None:
        raise ValueError("mamba executable not found")
    argv = [manager, "create", "--yes", "--prefix", str(prefix), "--file", str(lock)]
    print(json.dumps({"online_staging_argv": argv}), flush=True)
    # Conda prefix relocation prevents creating elsewhere then renaming. An
    # interrupted/failed installation remains for inspection; never call it ready.
    subprocess.run(argv, check=True, env=environment())
    check_records(prefix, packages(args.kind))
    print(f"Staged and package-verified: {prefix}")


def preflight(args):
    alignment, rprefix = args.alignment.resolve(), args.r_prefix.resolve()
    destination = new_destination(args.out)
    for kind, prefix in (("alignment", alignment), ("p0", rprefix)):
        check_lock(ROOT / f"config/{kind}-linux-64.explicit.txt", packages(kind))
        check_records(prefix, packages(kind))
    manifest = ROOT / "vendor/toolchain/manifest.json"
    source_info = json.loads(manifest.read_text())
    for item in source_info["retained_files"]:
        path = ROOT / item["path"]
        if digest(path) != item["sha256"]:
            raise ValueError(f"retained recipe/license hash mismatch: {path}")
    env = environment()
    env["PATH"] = str(alignment / "bin") + os.pathsep + env.get("PATH", "")
    versions = [("STAR", "--version", r"(?m)^2\.7\.10b\s*$"),
                ("hisat2", "--version", r"version 2\.2\.1(?:\s|$)"),
                ("hisat2-build", "--version", r"version 2\.2\.1(?:\s|$)"),
                ("samtools", "--version", r"(?m)^samtools 1\.17\s*$"),
                ("featureCounts", "-v", r"featureCounts v2\.0\.6(?:\s|$)")]
    # Only the requested new artifact directory is modified by preflight.
    scratch = Path(tempfile.mkdtemp(prefix=".preflight-", dir=destination.parent))
    try:
        observed = []
        for name, flag, pattern in versions:
            binary = alignment / "bin" / name
            result = subprocess.run([str(binary), flag], capture_output=True, text=True,
                                    errors="replace", timeout=30, env=env)
            output = result.stdout + result.stderr
            (scratch / f"{name}.version.txt").write_text(output)
            if result.returncode or not re.search(pattern, output):
                raise ValueError(f"{name}: pinned version probe failed: {output}")
            observed.append({"name": name, "executable": str(binary), "sha256": digest(binary),
                             "version_output": output})
        with (scratch / "R.log").open("w") as log:
            subprocess.run([str(rprefix / "bin/Rscript"), "--vanilla",
                            str(ROOT / "tools/preflight_r.R"), str(scratch)],
                           stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180, env=env)
        report = {"platform": platform.platform(), "alignment_prefix": str(alignment),
                  "r_prefix": str(rprefix), "tools": observed,
                  "tool_manifest_sha256": digest(manifest),
                  "locks": {kind: digest(ROOT / f"config/{kind}-linux-64.explicit.txt")
                            for kind in ("alignment", "p0")},
                  "scope": "Installed package records, versions, R computation and PNG; no installation or fetch."}
        (scratch / "verified.json").write_text(json.dumps(report, indent=2) + "\n")
        if destination.exists() or destination.is_symlink():
            raise ValueError("output appeared during preflight")
        scratch.rename(destination)
    except BaseException:
        print(f"Failed preflight diagnostics retained at {scratch}", file=sys.stderr)
        raise
    print(f"PASS: offline toolchain and R preflight: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    install = commands.add_parser("stage", help="explicit ONLINE installation; never invoked by compute")
    install.add_argument("kind", choices=("alignment", "p0"))
    install.add_argument("--prefix", type=Path, required=True)
    install.add_argument("--mamba", default="mamba")
    check = commands.add_parser("preflight", help="OFFLINE package/version/R/PNG checks; no install")
    check.add_argument("--alignment", type=Path, default=ROOT / ".deps/alignment")
    check.add_argument("--r-prefix", type=Path, default=ROOT / ".deps/p0")
    check.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        (stage if args.command == "stage" else preflight)(args)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
