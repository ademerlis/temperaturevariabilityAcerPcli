#!/bin/bash
#BSUB -J simplify_gtf
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o simplify_gtf.out
#BSUB -e simplify_gtf.err
#BSUB -u and128@miami.edu
#BSUB -N

# Specify the input and output files
input_gtf="path/to/your/input.gtf"
output_gtf="path/to/your/output.gtf"

# Use sed to simplify the identifiers
sed -e 's/gene_id "P5673_\(.*\)";/gene_id "\1";/g' \
    -e 's/transcript_id "P5673_\(.*\)";/transcript_id "\1";/g' \
    "$input_gtf" > "$output_gtf"
