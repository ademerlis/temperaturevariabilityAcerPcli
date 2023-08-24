#!/bin/bash
#BSUB -J fastqc
#BSUB -q bigmem
#BSUB -n 8
#BSUB -P and_transcriptomics
#BSUB -o fastqc_%J.out
#BSUB -e fastqc_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

cd ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/

fastqc *.fastq.gz --outdir ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/fastqc_files_rawsequences
