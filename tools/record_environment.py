#!/usr/bin/env python3
"""Record installed package provenance and post-link source archive hashes."""
import hashlib
import json
from pathlib import Path
import sys

prefix = Path(sys.argv[1]).resolve()
packages = []
for path in sorted((prefix / "conda-meta").glob("*.json")):
    record = json.loads(path.read_text())
    packages.append({key: record.get(key) for key in
                     ("name", "version", "build", "subdir", "url", "sha256", "md5", "license", "depends")})
archives = []
for path in sorted((prefix / "share").rglob("*.tar.gz*")):
    if path.is_file():
        archives.append({"path": str(path.relative_to(prefix)),
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
post_link_sources = {}
source_catalog = prefix / "share/bioconductor-data-packages/dataURLs.json"
if source_catalog.exists():
    catalog = json.loads(source_catalog.read_text())
    for record in packages:
        key = record["name"].removeprefix("bioconductor-") + "-" + record["version"]
        if key in catalog:
            post_link_sources[key] = catalog[key]
print(json.dumps({"platform": "linux-64", "packages": packages,
                  "post_link_archives": archives,
                  "post_link_sources": post_link_sources}, indent=2))
