# Retired workflow recovery

The current entry point is `bash run_pipeline.sh plan|local|submit CONFIG.tsv`.
It replaces the Python metadata/merge scripts, Snakemake DAG, generated YAML/Slurm
launcher, online R script and obsolete download/HTSeq/minimap shell pipelines.
Preparation and verification helpers now use R or shell.

For historical investigation, extract the pre-retirement tree into a separate,
empty directory (do not restore it over the current checkout):

```sh
mkdir -p /tmp/rnaseqdge-legacy
git archive 2f28776 | tar -x -C /tmp/rnaseqdge-legacy
```

That tree includes `metadata.py`, `merge_transcripts.py`, `pipeline.smk`,
`run_deg_analysis.R`, `requirements.txt`, old launchers and Python test helpers.
It requires its historical dependencies and is not the supported runtime.
Historical reports, locks and validation notes describe the software used at the
time; mentions of Python there do not declare a current compute dependency.

The retired helper-specific probes are `tools/p0_legacy_probe.py` and
`tests/integration/check_p0_legacy.py`. Their assertions concern bugs in the
retired merger, so they are preserved through the historical tree rather than
ported into the supported workflow. All remaining preparation/check helpers have
R or shell replacements with the same basename. `tools/toolchain.sh` additionally
provides a shell bootstrap before R is installed.
