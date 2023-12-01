#!/usr/bin/env bash
#BSUB -e bowtie2build_LocatelliShoguchiConcat.err
#BSUB -o bowtie2build_LocatelliShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8

workdir="/scratch/projects/and_transcriptomics"

bowtie2-build ${workdir}/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa, \
${workdir}/genomes/Symbiodinium/syma_transcriptome_37.fasta \
Locatelli_Shoguchi_concat
