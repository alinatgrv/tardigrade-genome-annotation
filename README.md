# Tardigrade Genome Annotation Pipeline

Chromatin-associated protein candidate identification in Ramazzottius varieornatus.

---

## Project overview

This repository contains a small reproducible annotation workflow combining:

- BLAST (homology search)
- HMMER (Pfam domain detection)
- Subcellular localization prediction
- Manual candidate prioritization

The goal is to demonstrate practical command-line bioinformatics and structured project organization.

---

## Repository structure

data/
  raw/                # original input files (not versioned)
  processed/          # cleaned inputs

results/
  blast/              # BLAST outputs
  hmmer/              # Pfam domain outputs
  localization/       # TargetP / WoLF PSORT outputs

report/
  figures/            # figures used in final report
  report.md           # narrative analysis

scripts/              # runnable bash scripts

---

## Environment

Recommended setup:

mamba create -n tardigrade -c conda-forge -c bioconda blast hmmer -y
conda activate tardigrade

---

## Example usage

bash scripts/run_blast.sh
bash scripts/run_hmmer.shgit ls-files | grep -E "\.pl$|\.pm$"

---

## Notes

This project is intended as a research/engineering demonstration of a small protein annotation workflow.
