#BSUB -u and128@miami.edu

#specify variables and paths

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads"


cd ${projdir}/take_4/trimmed_files

data=($(ls *.trim))

for samp in "${data[@]}" ; do \

#build script
echo "making fastqc script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_fastqc
#BSUB -e ${projdir}/logs/${samp}_fastqc.err
#BSUB -o ${projdir}/logs/${samp}_fastqc.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd ${projdir}/take_4/trimmed_files

module load fastqc/0.10.1
fastqc ${samp} --outdir ${projdir}/take_4/trimmed_files/fastqc_output/
echo \"Fastqc script of $samp submitted\"
" > ${projdir}/take_4/trimmed_files/${samp}_fastqc.job

bsub < ${projdir}/take_4/trimmed_files/${samp}_fastqc.job

done
