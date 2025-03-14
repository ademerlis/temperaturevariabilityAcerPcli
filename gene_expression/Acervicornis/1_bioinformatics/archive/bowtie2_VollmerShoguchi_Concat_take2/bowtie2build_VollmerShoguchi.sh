#!/usr/bin/env bash
#BSUB -e bowtie2build_VollmerShoguchiConcat.err
#BSUB -o bowtie2build_VollmerShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8

workdir="/scratch/projects/and_transcriptomics"

bowtie2-build ${workdir}/genomes/Acer_2023/Vollmer_2023/output.fasta, \
${workdir}/genomes/Symbiodinium/syma_transcriptome_37.fasta \
Vollmer_Shoguchi_concat
