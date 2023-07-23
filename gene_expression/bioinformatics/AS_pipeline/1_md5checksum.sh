#!/bin/bash
#BSUB -J transfer_checks
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o transfer_checks.out
#BSUB -e transfer_checks.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023"

# generate md5
md5sum ${and}/fastq_rawreads/*.gz > ${and}/tempvariability.md5

# count number of sequences per file

zcat ${and}/fastq_rawreads/*.fastq.gz | echo $((`wc -l`/4)) > ${and}/fastq_rawreads/rawread.counts.txt
