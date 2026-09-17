#!/usr/bin/env python3
"""Select every fifth record from verified 50k prefixes for the local deadline."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys

source, target = map(lambda x: Path(x).resolve(), sys.argv[1:])
target.mkdir()
(target / "reads").mkdir()
provenance = json.loads((source / "provenance.json").read_text())
for read in provenance["reads"]:
    original = source / read["path"]
    if hashlib.sha256(original.read_bytes()).hexdigest() != read["sha256"]:
        raise SystemExit("staged prefix checksum mismatch")
    count = 0
    with original.open("rb") as incoming, (target / read["path"]).open("wb") as outgoing:
        for index in range(50000):
            record = [incoming.readline() for _ in range(4)]
            if any(not line for line in record):
                raise SystemExit("short source prefix")
            if index % 5 == 0:
                outgoing.writelines(record)
                count += 1
        if incoming.read(1) or count != 10000:
            raise SystemExit("unexpected source/selected record count")
    read.update(source_prefix_records=50000, selected_records=count,
                prefix_records=count, selection="zero-based records 0,5,...,49995 of 50k archive prefix",
                source_prefix_sha256=read["sha256"],
                sha256=hashlib.sha256((target / read["path"]).read_bytes()).hexdigest(),
                bytes=(target / read["path"]).stat().st_size)
for name in ("samples.tsv", "runs.tsv", "analysis.tsv", "contrasts.tsv", "genes.tsv", "annotation.tsv"):
    shutil.copyfile(source / name, target / name)
with (source / "references.tsv").open() as incoming:
    rows = list(csv.DictReader(incoming, delimiter="\t"))
for row in rows:
    row["path"] = str(source / row["path"])
with (target / "references.tsv").open("w", newline="") as outgoing:
    writer = csv.DictWriter(outgoing, rows[0].keys(), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
provenance.update(records_per_run=10000, sampling="every fifth record of first50k archive prefix; nonrandom",
                  parent_stage=str(source), reason="50k STAR pilot exceeded budget for six workflows")
(target / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
(target / "STAGED").write_text("six verified 10k read selections for reduced local check\n")
print("STAGED six samples x 10,000 reads")
