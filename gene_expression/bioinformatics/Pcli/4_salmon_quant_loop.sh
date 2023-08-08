#BSUB -u and128@miami.edu

#specify variables and paths

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
echo "Aligning ${sample}"

echo '#!/bin/bash' > "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -q bigmem' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -J '"${sample}"_salmonquant'' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -o '"${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"$sample"_salmonquant%J.out'' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -e '"${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"$sample"_salmonquant%J.err'' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -n 8' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo '#BSUB -N' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job

echo ${and}/programs/salmon-1.5.2_linux_x86_64/bin/salmon quant -i \ ${and}/genomes/Pcli/Pcli_transcriptome_index -l U -r ${sample} \
--validateMappings \
-o ${and}/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/salmon_quant_files/"${sample}"_salmonquant >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job

echo 'echo' "Salmon quantification of $sample complete"'' >> "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job
echo "Salmon quantification script of $sample submitted"
bsub < "${and}"/Ch2_temperaturevariability2023/AS_pipeline/4_Pcli_specific/scripts/"${sample}"_salmonquant.job ; \
done
