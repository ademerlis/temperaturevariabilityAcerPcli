#!/bin/bash
#BSUB -J trim_stresshardening
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_stresshardeningA%J.out
#BSUB -e trim_stresshardeningA%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

conda activate condaenv

and="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq"

for sample in ${and}/*.gz ; 

do \

tagseq_clipper.pl $sample | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${sample/.gz}.trim ; \

done
