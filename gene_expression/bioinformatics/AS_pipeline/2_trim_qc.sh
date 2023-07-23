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

for sample in ${and}/Ch2_temperaturevariability2023/fastq_rawreads/*.fastq.gz ;
do \
${and}/programs/fastp --in1 ${sample} --out1 ${sample}.clean.processed --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 ; \
done

cd ${and}/Ch2_temperaturevariability2023/fastq_rawreads/

fastqc *.clean.processed

multiqc .
