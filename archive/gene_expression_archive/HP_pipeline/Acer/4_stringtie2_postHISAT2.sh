#!/bin/bash
#BSUB -J stringtie2_postHISAT2
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie2_postHISAT2%J.out
#BSUB -e stringtie2_postHISAT2%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Acer_aligned_bam_files"

data=($(ls *.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -p 8 -e -B -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
echo "StringTie assembly for seq file ${i}" $(date) ; \
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date) 
