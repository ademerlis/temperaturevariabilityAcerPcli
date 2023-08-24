#!/bin/bash
#BSUB -J Pcli_transcriptome_salmon_index
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o Pcli_transcriptome_salmon_index%J.out
#BSUB -e Pcli_transcriptome_salmon_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon index -t ${and}/genomes/Pcli/clean_Pcli_transcriptome_final.fasta -i ${and}/genomes/Pcli/Pcli_transcriptome_index
