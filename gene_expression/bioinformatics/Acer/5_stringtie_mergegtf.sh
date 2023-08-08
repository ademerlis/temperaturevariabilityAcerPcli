#!/bin/bash
#BSUB -J stringtie_mergegtf
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_mergegtf%J.out
#BSUB -e stringtie_mergegtf%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Acer_aligned_bam_files"

ls *.gtf > gtf_list.txt

${and}/programs/stringtie-2.2.1/stringtie --merge -e -p 8 -G ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -o stringtie_acerv_merged.gtf gtf_list.txt
echo "Stringtie merge complete" $(date)

/scratch/projects/and_transcriptomics/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare -r ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -G -o merged stringtie_acerv_merged.gtf
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python ${and}/programs/stringtie-2.2.1/prepDE.py -g gene_count_acerv_matrix.csv -i listGTF.txt 

echo "Gene count matrix compiled." $(date)
