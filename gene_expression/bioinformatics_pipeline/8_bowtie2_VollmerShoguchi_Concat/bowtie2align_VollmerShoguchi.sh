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
#BSUB -J ${samp}_bowtie2align_VollmerShoguchi
#BSUB -e ${projdir}/bowtie2_VollmerShoguchi_concat/logs/${samp}_bowtie2align_VollmerShoguchi.err
#BSUB -o ${projdir}/bowtie2_VollmerShoguchi_concat/logs/${samp}_bowtie2align_VollmerShoguchi.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer\"

bowtie2 --local -U ${samp} -x Vollmer_Shoguchi_concat --un ${samp}.unaligned -k 5 -S ${samp}.sam

done

" > ${projdir}/bowtie2_VollmerShoguchi_concat/${samp}_bowtie2align_VollmerShoguchi.job

bsub < ${projdir}/bowtie2_VollmerShoguchi_concat/${samp}_bowtie2align_VollmerShoguchi.job

done
