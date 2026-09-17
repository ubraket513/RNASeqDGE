#!/usr/bin/env python3
"""Read-only local tool inventory; never installs tools or submits jobs."""
import importlib.util
import json
import platform
import shutil
import subprocess
import sys


def probe(argv):
    try:
        result = subprocess.run(argv, capture_output=True, text=True, timeout=30)
        return {"exit_code": result.returncode,
                "output": (result.stdout + result.stderr).strip()}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"error": str(exc)}


def main():
    flags = {"g++": "--version", "make": "--version", "Rscript": "--version",
             "STAR": "--version", "hisat2": "--version", "featureCounts": "-v",
             "samtools": "--version", "snakemake": "--version",
             "fasterq-dump": "--version", "sbatch": "--version", "node": "--version"}
    report = {"platform": platform.platform(), "python": sys.version,
              "python_executable": sys.executable, "tools": {}}
    for name, flag in flags.items():
        path = shutil.which(name)
        report["tools"][name] = {"path": path, **(probe([path, flag]) if path else {})}
    report["python_modules"] = {
        name: importlib.util.find_spec(name) is not None
        for name in ("pandas", "biomart", "loguru", "Bio", "snakemake")}
    if shutil.which("Rscript"):
        report["r_packages"] = probe([
            "Rscript", "--vanilla", "-e",
            'for (p in c("DESeq2","apeglm","GEOquery","optparse","tidyverse",'
            '"EnhancedVolcano","RColorBrewer","pheatmap","ggrepel",'
            '"org.Hs.eg.db","BiocParallel")) cat(p, '
            'if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) '
            'else "MISSING", "\\n")'])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
