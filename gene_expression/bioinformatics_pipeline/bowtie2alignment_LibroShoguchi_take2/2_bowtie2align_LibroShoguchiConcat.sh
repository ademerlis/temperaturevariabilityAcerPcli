#!/usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"


#BSUB -e bowtie2align_LibroShoguchiConcat.err
#BSUB -o bowtie2align_LibroShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -f -U "$f" -x Host_concat --keep_unal -k 5 -S "$output_file"
done


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -U "$f" -x Host_concat --un "$output_file".unaligned -k 5 -S "$output_file"
done
