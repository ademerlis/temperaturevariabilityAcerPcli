#!/bin/bash
#BSUB -J countreads_trim
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o countreads_trim.out
#BSUB -e countreads_trim.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in trimmed fastq files

#specify variables and paths

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4"

output_file="countreads_results.txt"

# Default file pattern
glob="\.trim"

# Check if an argument is provided
if [ "$1" ]; then
    glob="$1"
fi

# Loop through files matching the pattern
for f in *$glob*; do
    # Count the number of lines in the file
    nrd=$(cat "$f" | wc -l)

    # Divide the line count by 4
    nrd=$((nrd / 4))

    # Print the filename and the calculated number
    echo -e "$f\t$nrd" >> "$output_file"
done

echo "Results have been saved to $output_file"
