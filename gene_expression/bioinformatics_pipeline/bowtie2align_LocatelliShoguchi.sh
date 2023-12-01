#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

data=($(ls *.trim))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_bowtie2align_LocatelliShoguchi
#BSUB -e ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.err
#BSUB -o ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer\"

bowtie2 --local -U ${samp} -x Locatelli_Shoguchi_concat --un ${samp}.unaligned -k 5 -S ${samp}.sam

done

" > ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

bsub < ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

done
