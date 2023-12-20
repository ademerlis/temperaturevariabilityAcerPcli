#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/symbionts"

data=($(ls *.sym.clean))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_align_sym_final
#BSUB -e ${projdir}/logs/${samp}_align_sym_final.err
#BSUB -o ${projdir}/logs/${samp}_align_sym_final.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/symbionts\"

bowtie2 --local -U ${samp} -x ../Syma_index -S ${samp}.sym.clean.sam --no-hd --no-sq --no-unal

" > ${projdir}/${samp}_align_sym_final.job

bsub < ${projdir}/${samp}_align_sym_final.job

done
