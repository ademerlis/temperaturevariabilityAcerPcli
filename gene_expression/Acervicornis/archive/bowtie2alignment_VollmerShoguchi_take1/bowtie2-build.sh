#!/bin/bash
#BSUB -J bowtie2-build
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o bowtie2-build.out
#BSUB -e bowtie2-build.err
#BSUB -u and128@miami.edu
#BSUB -N

workdir="/scratch/projects/and_transcriptomics/genomes"

bowtie2-build ${workdir}/Acer_2023/GCA_032359415.1_NEU_Acer_K2_genomic.fna, \
${workdir}/Symbiodinium/syma_transcriptome_37.fasta,\
${workdir}/Breviolum/Symb_Dip.fna,\
${workdir}/Cladocopium/C124.annotated.fa,\
${workdir}/Durusdinium/102_symbd_transcriptome_nucl.fa \
Host_concat
