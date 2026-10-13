#!/bin/sh
set -eu
# Snakemake runs this from the project root after creating a rule environment.
python scripts/install_vep.py
