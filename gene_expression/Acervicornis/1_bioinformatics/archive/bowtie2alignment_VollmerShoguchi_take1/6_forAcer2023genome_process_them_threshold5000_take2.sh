#!/usr/bin/env bash
#BSUB -P and_transcriptomics
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.err
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/UTR_add.out

module load python/3.8.7
cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/
GTF_FILE="Acer_parsed.gtf"

mkdir temp_files;
mkdir results;
grep '#' ${GTF_FILE} > ${GTF_FILE}_hashtag_lines
# split input gtf into fw and rev
awk -F'\t' '$7=="+" {print $0}' ${GTF_FILE} > ${GTF_FILE}_fw.gtf
awk -F'\t' '$7=="-" {print $0}' ${GTF_FILE} > ${GTF_FILE}_rev.gtf

cd /scratch/projects/and_transcriptomics/genomes/Acer_2023/UTR_add_extend_GTF-main/ ;
module load python/3.8.7 ;
python3 Extend_by_threshold_fw.py --input ${GTF_FILE}_fw.gtf --threshold 5000 > temp_files/${GTF_FILE}_fw.gtf ;
python3 Extend_by_threshold_rev.py --input ${GTF_FILE}_rev.gtf --threshold 5000 > temp_files/${GTF_FILE}_rev.gtf ;
cat temp_files/${GTF_FILE}_fw.gtf > temp_files/unsorted ;
cat temp_files/${GTF_FILE}_rev.gtf >> temp_files/unsorted ;
sort -k1,1V -k4,4n -k5,5rn temp_files/unsorted > temp_files/sorted ;
cat ${GTF_FILE}_hashtag_lines > results/${GTF_FILE}_ext_5000.gtf ;
cat temp_files/sorted >> results/${GTF_FILE}_ext_5000.gtf
