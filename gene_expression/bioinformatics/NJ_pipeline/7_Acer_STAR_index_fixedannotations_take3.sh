#!/bin/bash
#BSUB -J Acer_star_index_fixedannotations_take3
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take3 \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 13
