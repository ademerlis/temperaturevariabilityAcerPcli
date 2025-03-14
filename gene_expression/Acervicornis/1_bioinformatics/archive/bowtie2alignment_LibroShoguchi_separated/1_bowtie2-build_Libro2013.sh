#!/bin/bash
#BSUB -J bowtie2-build_Libro2013
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o bowtie2-build_Libro2013.out
#BSUB -e bowtie2-build_Libro2013.err
#BSUB -u and128@miami.edu
#BSUB -N

workdir="/scratch/projects/and_transcriptomics/genomes"

bowtie2-build ${workdir}/Acer_2023/Libro_2013/acer.fasta \
Acer_index
