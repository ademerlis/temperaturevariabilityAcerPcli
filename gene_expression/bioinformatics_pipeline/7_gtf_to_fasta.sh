#!/usr/bin/env bash

# convert gtf to fasta
module load cufflinks/2.2.1
gffread -w output.fasta -g /scratch/projects/and_transcriptomics/genomes/Acer_2023/Vollmer_2023/GCA_032359415.1_NEU_Acer_K2_genomic.fna /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/results/Acer_parsed.gtf_ext_5000.gtf
