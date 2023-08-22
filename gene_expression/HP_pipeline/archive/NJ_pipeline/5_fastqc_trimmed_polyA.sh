#!/bin/bash
#~/scripts/fastqc_trimmed.sh
#purpose: quality checking of trimmed reads using FASTQC on Pegasus compute node

#BSUB -J fastqc_trimmed
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_trimmed.out
#BSUB -e fastqc_trimmed.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq"

module load fastqc/0.10.1

cd ${and}/trimmed

fastqc *.fastq.gz
--outdir ${and}/trimmed/


