#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

data=($(ls *.host))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_align_Host
#BSUB -e ${projdir}/logs/${samp}_align_Host.err
#BSUB -o ${projdir}/logs/${samp}_align_Host.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer\"

bowtie2 --local -U ${samp} -x Acer_index -S ${samp}.host.sam --no-hd --no-sq --no-unal --al ./${samp}.fastq.host.clean --un junk/${samp}.fastq.host.un

" > ${projdir}/${samp}_align_Host.job

bsub < ${projdir}/${samp}_align_Host.job

done
