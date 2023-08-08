#!/bin/bash
#BSUB -J HISAT2
#BSUB -q bigmem
#BSUB -n 16
#BSUB -P and_transcriptomics
#BSUB -o HISAT2%J.out
#BSUB -e HISAT2%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load samtools/1.3
module load python/3.8.7

/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2-build -f ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta ${and}/genomes/Acer/Acer_reference_genome_hisat2
echo "Reference genome indexed. Starting alignment" $(date)

cd /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/
array=($(ls *.fastq.gz))
for i in ${array[@]};
 do \
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	/scratch/projects/and_transcriptomics/programs/hisat2-2.2.1/hisat2 -p 8 --dta -x ${and}/genomes/Acer/Acer_reference_genome_hisat2 -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam ; 
done
