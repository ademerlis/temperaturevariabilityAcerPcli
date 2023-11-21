#!/usr/bin/env bash
#BSUB -e bowtie2map.err
#BSUB -o bowtie2map.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8


name_position=$1

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

glob="\.trim"

if [ ! -z "$1" ]; then
    glob="$1"
fi
# Loop through files matching the pattern
for f in *$glob*; do
    if [ -n "$name_position" ]; then
        # Split the filename and extract the part based on the provided position
        IFS='._' read -ra parts <<< "$f"
        outname="${parts[$((name_position - 1))]}.sam"
    else
        # Use the entire filename if no position is provided
        outname="${f}.sam"
    fi
    bowtie2 --local -x Host_concat -U $f -S $outname --no-hd --no-sq --no-unal -k 5
done
