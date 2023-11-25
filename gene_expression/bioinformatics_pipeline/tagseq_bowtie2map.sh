#!/usr/bin/env bash
#BSUB -e bowtie2map_libro.err
#BSUB -o bowtie2map_libro.out
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
