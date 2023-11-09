#!/bin/bash
#BSUB -J trim_fastp.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_fastp_%J.out
#BSUB -e trim_fastp_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.gz))

for sample in ${data[@]} ;
do \
${and}/programs/fastp --in1 ${sample} --out1 ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -h report_${sample}.html --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
fastqc ${and}/Ch2_temperaturevariability2023/fastp_processed/clean.${sample} -o ${and}/Ch2_temperaturevariability2023/trimmed_qc_files/ ; \
done

echo "Read trimming of adapters complete."

multiqc ${and}/Ch2_temperaturevariability2023/trimmed_qc_files/

mv multiqc_report.html trimmed_qc_files/
mv multiqc_data trimmed_qc_files/

echo "Cleaned MultiQC report generated."
