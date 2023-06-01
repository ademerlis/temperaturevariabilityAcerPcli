#!/bin/bash
#~/scripts/fastqc_stresshardening2022.sh
#/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/scripts/fastqc_stresshardening2022.sh
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J stresshardening2022_fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_stresshardening2022.out
#BSUB -e fastqc_stresshardening2022.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

module load fastqc/0.10.1

and="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq"

cd ${and}
fastqc *.fastq.gz
--outdir ${and}/fastqc/
