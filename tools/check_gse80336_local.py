#!/usr/bin/env python3
"""Run the explicitly reduced, serial real-read check under a wall-time deadline."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seconds", type=int, default=1800)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--index-cache", type=Path)
    args = parser.parse_args()
    if not 1 <= args.seconds <= 1800 or not 1 <= args.replicates <= 3:
        parser.error("seconds must be 1..1800 and replicates 1..3")
    args.inputs = args.inputs.resolve()
    args.out = args.out.absolute()
    args.out.mkdir()  # Refuse overwrite; preserve incomplete checks for diagnosis.
    deadline = time.monotonic() + args.seconds
    native = ROOT / "build/rnaseq"
    native_hash = hashlib.sha256(native.read_bytes()).hexdigest()
    records = []
    def record_progress():
        (args.out / "runs.json").write_text(json.dumps(dict(
            scope="6 samples, sampled archive-prefix reads, chromosome22 only; not whole-genome inference",
            native_sha256=native_hash, threads=2, workers=1,
            requested_replicates=args.replicates,
            complete=len(records) == 2 * args.replicates and all(r["status"] == 0 for r in records),
            runs=records), indent=2) + "\n")
    record_progress()
    env = dict(os.environ, PATH="/usr/bin:/bin", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", R_LIBS_USER="/dev/null", R_LIBS_SITE="/dev/null")
    for key in ("R_HOME", "R_LIBS", "LD_LIBRARY_PATH", "LD_PRELOAD"):
        env.pop(key, None)
    for repeat in range(1, args.replicates + 1):
        # Alternate order to reduce a consistent first-backend cache advantage.
        backends = ("hisat2", "star") if repeat % 2 else ("star", "hisat2")
        for backend in backends:
            remaining = int(deadline - time.monotonic())
            if remaining <= 0:
                raise SystemExit("local-check deadline reached; partial results retained")
            if hashlib.sha256(native.read_bytes()).hexdigest() != native_hash:
                raise SystemExit("native executable changed during check")
            name = f"{backend}-{repeat}"
            settings = {k: args.inputs / f"{k}.tsv" for k in
                        ("samples", "runs", "references", "analysis", "contrasts", "genes")}
            settings.update(bin_dir=ROOT / ".deps/alignment/bin", rscript=ROOT / ".deps/p0/bin/Rscript",
                            r_script=ROOT / "run_deg_analysis_offline.R",
                            tool_lock=ROOT / "config/alignment-linux-64.explicit.txt",
                            r_lock=ROOT / "config/p0-linux-64.explicit.txt", run_dir=args.out / name,
                            index_cache=args.index_cache.resolve() if args.index_cache else args.out / "index-cache", backend=backend, threads=2, workers=1,
                            star_sa_bases=11, star_chr_bits=16)
            config = args.out / f"{name}.tsv"
            with config.open("w", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(["key", "value"])
                writer.writerows(settings.items())
            command = ["/usr/bin/timeout", "--signal=TERM", "--kill-after=15s", f"{remaining}s",
                       "/usr/bin/time", "-v", "-o", str(args.out / f"{name}.time.txt"),
                       str(native), "workflow-local", str(config)]
            print(f"START {name}: {remaining}s remain", flush=True)
            start = time.monotonic()
            with (args.out / f"{name}.log").open("w") as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env,
                                           start_new_session=True)
                try:
                    status = process.wait()
                except BaseException:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    raise
            record = dict(backend=backend, repeat=repeat, status=status,
                          wall_seconds=time.monotonic() - start, command=command,
                          index_condition="reuse from pilot" if args.index_cache else ("build" if repeat == 1 else "reuse"),
                          cache_condition="uncontrolled OS cache; no cache eviction requested")
            records.append(record)
            record_progress()
            print(f"END {name}: exit={status}, {record['wall_seconds']:.1f}s", flush=True)
            if status != 0:
                raise SystemExit(f"{name} failed; see retained log and stage diagnostics")
    print("All requested reduced workflows completed; scientific comparison remains required.")


if __name__ == "__main__":
    main()
