# Juicer course version

This repository contains a pinned Juicer copy for the
`OMICS_course_spring_2026` Hi-C practical class.

It is based on Aiden Lab Juicer:
https://github.com/aidenlab/juicer

The goal of this repository is reproducibility for the course: students clone
the same tested Juicer layout and run the same commands during the practical.
The original license is kept in `LICENSE`.

## What is included

- `scripts/juicer.sh`
- helper scripts from `scripts/common/`
- `juicer_tools.jar` used by this Juicer layout
- helper scripts from `misc/`, including `generate_site_positions.py`

Large reference files and precomputed restriction-site files are not included.
For the course, restriction sites are generated from the selected reference
genome with `misc/generate_site_positions.py`.

## Installation

From the course practice directory:

```bash
git clone https://github.com/dpanc2/juicer_course_version.git tools/juicer
```

## Example run

```bash
bash tools/juicer/scripts/juicer.sh \
  -D "$(pwd)/tools/juicer" \
  -d "$(pwd)/data/juicer/MoPh7" \
  -g T2T_human \
  -z "$(pwd)/data/reference/T2T_human.fna" \
  -p "$(pwd)/data/reference/chrom.sizes" \
  -y "$(pwd)/data/reference/restriction_sites_DpnII.txt" \
  -s DpnII \
  -t 4
```

The expected output is:

```text
data/juicer/MoPh7/aligned/inter_30.hic
```
