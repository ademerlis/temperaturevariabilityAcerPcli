#!/usr/bin/env bash
#BSUB -e bowtie2map_libro_3.err
#BSUB -o bowtie2map_libro_3.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob=".trim"

# Loop through files matching the pattern
for f in *$glob*; do
    output_file="${f%$glob}.sam"
    bowtie2 --local -U "$f" -x Host_concat --un "$output_file".unaligned -k 5 -S "$output_file"
done
