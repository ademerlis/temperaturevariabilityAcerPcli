#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_quant
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -o Pcli_transcriptome_salmon_quant%J.out
#BSUB -e Pcli_transcriptome_salmon_quant%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${and}/genomes/Pcli/Pcli_transcriptome_index -l U -r ${sample} \
--validateMappings --gcBias --reduceGCMemory -o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files ; \
done
