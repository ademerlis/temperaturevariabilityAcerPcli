#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.fastq.gz))

for samp in "${data[@]}" ; do \

#build script
echo "making cutadapt script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_trim
#BSUB -e ${projdir}/logs/${samp}_trim.err
#BSUB -o ${projdir}/logs/${samp}_trim.out
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads\"

echo \"Gunzipping file  ${samp} then removing adaptors and low quality reads...\"

gunzip -c ${samp} > ${samp}_temp.fastq
${and}/Ch2_temperaturevariability2023/0_scripts/tagseq_clipper.pl ${samp}_temp.fastq | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${projdir}/${samp/.fastq.gz/}.trim

" > ${projdir}/${samp}_trim.job

#submit script

bsub < ${projdir}/${samp}_trim.job

done

echo "All done!"
