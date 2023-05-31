#!/bin/bash
#~/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/scripts/fastqc/fastqc_stresshardening2022.job
#/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/scripts/fastqc/fastqc_stresshardening2022.job
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_stresshardening2022.out
#BSUB -e fastqc_stresshardening2022.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/" 

cd ${and}
for SAMP in *.fastq.gz

do

module load java/1.8.0_60
module load  \
${and}/$SAMP \
--outdir ${and}/fastqc_results
done
