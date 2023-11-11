#!/bin/bash
#BSUB -J trim_test_2
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o trim_test_2.out
#BSUB -e trim_test_2.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

${and}/programs/TrimGalore-0.6.10/trim_galore temp.fastq
--adapter2 "AAAAAAAA" \
--adapter "AGATCGG" \
--quality 15 \
--length 25 \ 
-o ${and}.temp.fastq.trim
