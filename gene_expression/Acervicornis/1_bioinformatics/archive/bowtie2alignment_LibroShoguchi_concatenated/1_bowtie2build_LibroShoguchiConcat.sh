#!/usr/bin/env bash
#BSUB -e bowtie2map_LibroShoguchiConcat.err
#BSUB -o bowtie2map_LibroShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8

workdir="/scratch/projects/and_transcriptomics/genomes"

bowtie2-build ${workdir}/Acer_2023/Libro_2013/acer.fasta, \
${workdir}/Symbiodinium/syma_transcriptome_37.fasta \
Host_concat
