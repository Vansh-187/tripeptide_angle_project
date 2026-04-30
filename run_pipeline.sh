#!/bin/bash

CORES=8


export PDB_DIR="${1:-pdbs}"
export STRIDE_BIN="${2:-/mnt/c/Users/VANSH/stride/stride}"


echo " PDB to STRIDE"
snakemake -s secondary_structure_analysis.smk --cores $CORES --keep-going

echo "STRIDE to Contexts"
snakemake -s create_context_tsv.smk --cores $CORES --keep-going

echo "Contexts to Angle list"
snakemake -s angle_pipeline.smk --cores $CORES --keep-going

echo "Plot"
Rscript scripts/plot_angles.R

echo "Done. Output: final/angle_plot.png"
