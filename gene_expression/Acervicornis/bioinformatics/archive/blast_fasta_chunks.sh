#BSUB -u and128@miami.edu

#specify variables and paths

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/genomes/Acer/annotatingTranscriptomes-master"

cd ${projdir}

data=($(ls subset*))

for file in "${data[@]}" ; do \

#build script
echo "making blast uniprot script for ${file}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${file}_blast_uniprot
#BSUB -e ${projdir}/logs/${file}_blast_uniprot.err
#BSUB -o ${projdir}/logs/${file}_blast_uniprot.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd ${projdir}

module load blast/2.2.29+

blastx -query ${file} -db uniprot_sprot.fasta -evalue 0.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out ${file}.br" > ${projdir}/${file}_blast_uniprot.job

bsub < ${projdir}/${file}_blast_uniprot.job

done
