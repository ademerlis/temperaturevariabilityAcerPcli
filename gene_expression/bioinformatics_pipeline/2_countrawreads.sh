#!/bin/bash
#BSUB -J countrawreads
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o countrawreads.out
#BSUB -e countrawreads.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in a bunch of fastq files

#specify variables and paths

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

output_file="countreads_results.txt"

glob=".fastq.gz"
if [ ! -z "$1" ]; then
    glob="$1"
fi

fqs=(*$glob)
for f in "${fqs[@]}"; do
    gunzip -c "$f" > "temp.fastq"  # Decompress the file to a temporary file
    nrd=$(cat "temp.fastq" | wc -l)
    nrd=$((nrd / 4))
    echo "$f	$nrd"
    echo "$f	$nrd" >> "$output_file"  # Append the results to the output file
    rm "temp.fastq"  # Remove the temporary file
done

echo "Results have been saved to $output_file"
