#!/bin/bash
#BSUB -J trim_qc.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_qc_%J.out
#BSUB -e trim_qc_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

cd ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/
array1=($(ls *.fastq.gz)) 
for sample in ${array1[@]} ;
do \
${and}/programs/fastp --in1 ${sample} --out1 clean.${sample} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 \
fastqc clean.${i} ; \
done

cd ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/

multiqc .
