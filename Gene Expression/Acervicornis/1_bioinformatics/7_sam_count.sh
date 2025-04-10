#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir=

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

bsub < ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

done
