#!/bin/bash
#BSUB -J bowtie2-build_Shoguchi2018
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o bowtie2-build_Shoguchi2018.out
#BSUB -e bowtie2-build_Shoguchi2018.err
#BSUB -u and128@miami.edu
#BSUB -N

workdir="/scratch/projects/and_transcriptomics/genomes"

bowtie2-build ${workdir}/Symbiodinium/syma_transcriptome_37.fasta \
Syma_index
