#!/bin/bash
#BSUB -J stringtie_genecounts
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_genecounts%J.out
#BSUB -e stringtie_genecounts%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/genecounts"

F="/scratch/projects/and_transcriptomics/Allyson_stresshardening_RNAseq/aligned/stringtie_reestimatedmerged_gtf_files"

data=($(ls *merge.gtf))

for i in ${data[@]} ;

do \
echo "${i} $F${i}" >> sample_list_acerv.txt ; \
done

python ${and}/programs/stringtie-2.2.1/prepDE.py -g gene_count_acerv_matrix.csv -i sample_list_acerv.txt
