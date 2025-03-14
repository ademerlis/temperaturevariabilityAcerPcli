#BSUB -u and128@miami.edu

and="/scratch/projects/and_transcriptomics"
cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.gz))

for sample in ${data[@]} ;
do \
echo '#!/bin/bash' > "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -q bigmem' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -J '"${sample}"_fastp'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -o '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"$sample"_fastp%J.out'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -e '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"$sample"_fastp%J.err'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -n 8' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '#BSUB -N' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo 'module load fastqc/0.10.1' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo '/scratch/projects/and_transcriptomics/programs/fastp '--in1 "${sample}" --out1 "${and}"/Ch2_temperaturevariability2023/fastp_processed/clean."${sample}" -h report_"${sample}".html -j "${sample}"_fastp.json --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
echo 'fastqc '"${and}"/Ch2_temperaturevariability2023/fastp_processed/clean."${sample}" -o "${and}"/Ch2_temperaturevariability2023/trimmed_qc_files/'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job
bsub < "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/"${sample}"_fastp.job ; \
done
