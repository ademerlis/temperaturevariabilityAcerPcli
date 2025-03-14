#!/bin/bash
#BSUB -J trim_SH
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -o trim_SH.out
#BSUB -e trim_SH.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_stresshardening_RNAseq/trimmed/trimmed_fastq_files/*.gz ;

do \
${and}/programs/TrimGalore-0.6.10/trim_galore ${sample} \
--polyA ; \

done

cd ${and}/Allyson_stresshardening_RNAseq/scripts

for f in *.fq.gz;
do \
mv -v -- "$f" "${f%.fq.gz}.fastq.gz"; \

done

module load fastqc/0.10.1

fastqc *.fastq.gz

multiqc *
