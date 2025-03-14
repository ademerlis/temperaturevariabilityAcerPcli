#BSUB -u and128@miami.edu

#specify variables and paths

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
echo "FastQC ${sample}"

echo '#!/bin/bash' > "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -q bigmem' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -J '"${sample}"_fastqc'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -o '"${and}"/Ch2_temperaturevariability2023/AS_pipeline/2_trimmed_reads/Acer_fastq_files/"$sample"_fastqc%J.out'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -e '"${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"$sample"_fastqc%J.err'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -n 8' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo '#BSUB -N' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo 'module load fastqc/0.10.1
fastqc '"${sample}" --outdir /scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo 'echo' "Fastqc of $sample complete"'' >> "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job
echo "Fastqc script of $sample submitted"
bsub < "${and}"/Ch2_temperaturevariability2023/2_trimmed_reads/Acer_fastq_files/"${sample}"_fastqc.job ; \
done
