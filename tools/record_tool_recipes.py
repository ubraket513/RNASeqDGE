#!/usr/bin/env python3
"""Retain pinned conda source/build/license evidence (staging only; needs PyYAML)."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import yaml

ROOT = Path(__file__).resolve().parents[1]


def main():
    prefix = Path(sys.argv[1]).resolve()
    target = ROOT / "vendor/toolchain"
    if target.exists():
        raise SystemExit(f"refusing existing output: {target}")
    records = [json.loads(p.read_text()) for p in sorted((prefix / "conda-meta").glob("*.json"))]
    manifest = {"format": 1, "platform": "linux-64",
                "distribution": "Pinned conda binaries; recipes retained as build provenance, not rebuilt locally.",
                "transitive_dependencies": "config/alignment-linux-64.provenance.json",
                "source_hash_scope": "Upstream archive hashes declared by retained recipes; upstream archives not downloaded by this recorder.",
                "packages": [], "retained_files": []}
    for name in ("star", "hisat2", "subread", "samtools", "htslib"):
        record = next(r for r in records if r["name"] == name)
        cache = Path(record["extracted_package_dir"])
        archive = Path(record["package_tarball_full_path"])
        with archive.open("rb") as stream:
            checksum = hashlib.file_digest(stream, "sha256").hexdigest()
        if checksum != record["sha256"]:
            raise ValueError(f"package archive hash mismatch: {archive}")
        info = cache / "info"
        recipe = yaml.safe_load((info / "recipe/meta.yaml").read_text())
        dest = target / name
        dest.mkdir(parents=True)
        # Preserve recipe patch files, build environment pins and all license texts.
        shutil.copytree(info / "recipe", dest / "recipe")
        if (info / "licenses").is_dir():
            shutil.copytree(info / "licenses", dest / "licenses")
        for filename in ("about.json", "hash_input.json", "index.json", "git"):
            if (info / filename).is_file():
                shutil.copyfile(info / filename, dest / filename)
        manifest["packages"].append({
            **{key: record[key] for key in ("name", "version", "build", "url", "sha256", "license", "depends")},
            "source": recipe["source"],
            "recipe": str((dest / "recipe/meta.yaml").relative_to(ROOT)),
            "build_script": str((dest / "recipe/build.sh").relative_to(ROOT)),
            "build_flags_scope": "Literal upstream recipe scripts and conda_build_config.yaml retained; no local rebuild claim."})
    for path in sorted(target.rglob("*")):
        if path.is_file():
            manifest["retained_files"].append({"path": str(path.relative_to(ROOT)),
                                              "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Recorded {len(manifest['packages'])} verified package archives and {len(manifest['retained_files'])} provenance files")


if __name__ == "__main__":
    main()
