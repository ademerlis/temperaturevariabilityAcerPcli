# Bioinformatics pipelines for *A.cervicornis* and *P.clivosa* samples

Script written by: DeMerlis 

Last updated: 20230807

For each step of this analysis, I followed the pipelines of [Jill Ashey](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md?plain=1), [Dr. Ariana Huffmyer](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md), [Dr. Sam Gurr](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md), and [Zoe Dellaert](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/ZD_Heron-Pdam-gene-expression.md). All code was adapted from them. 

## 1) Download Data

I downloaded the BaseSpace GUI and downloaded the .fastq files for Acer and Pcli onto an external hard drive. Then, I uploaded them to the UM HPC Pegasus using 'scp'. 

## 2) QC Raw Reads

Can run FastQC module from Pegasus (version 0.10.1). 

```{bash}
#!/bin/bash
#BSUB -J fastqc
#BSUB -q bigmem
#BSUB -n 8
#BSUB -P and_transcriptomics
#BSUB -o fastqc_%J.out
#BSUB -e fastqc_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

cd ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads

fastqc *.fastq.gz --outdir ${and}/Ch2_temperaturevariability2023/2_fastqc_files_rawsequences
```

Install multiqc locally:
```{bash}
module load py-pip/20.2
pip install multiqc

nano ~/.bash_profile
export PATH=$PATH:/nethome/and128/.local/bin

#then save it and exit nano

source .bash_profile
#this runs the bash profile to update the paths
```

Then cd into 2_fastqc_files_rawsequences and run `multiqc .`

I will transfer the multiqc_report.html to my local drive so I can open it and view it.

this code worked from transferring from pegasus to local (first navigated on local to the folder i wanted the report to go in):
```{bash}
Allysons-MacBook-Pro-2: allysondemerlis$ scp and128@pegasus.ccs.miami.edu:/scratch/projects/and_transcriptomics/multiqc_report.html .
```

## 3) Trim Files Using Fastp

First, install Fastp locally:
```{bash}
cd programs
# download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

Then, run script to trim fastq.gz files and rename to "clean.{sample}.fastq.gz". 

"This script uses fastp to: 
- remove TagSeq-specific adapter sequences (--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)
- enable polyX trimming on 3' end at length of 6 (--trim_poly_x 6)
- filter by minimum phred quality score of >30  (-q 30)
- enable low complexity filter (-y)
- set complexity filter threshold of 50% required (-Y 50)" (explained by Dr. Sam Gurr)

```{bash}
#!/bin/bash
#BSUB -J trim.sh
#BSUB -q bigmem
#BSUB -n 8
#BSUB -R "rusage[mem=10000]"
#BSUB -P and_transcriptomics
#BSUB -o trim_%J.out
#BSUB -e trim_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd ${and}/Ch2_temperaturevariability2023/AS_pipeline/1_fastq_rawreads/
array1=($(ls *.fastq.gz)) 
for sample in ${array1[@]} ;
do \
${and}/programs/fastp --in1 ${sample} --out1 clean.${sample} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 ; \
done
```

## 4) QC Trimmed Reads

```{bash}
#!/bin/bash
#BSUB -J fastqc_Pclitrimmed
#BSUB -q bigmem
#BSUB -n 8
#BSUB -P and_transcriptomics
#BSUB -o fastqc_%J.out
#BSUB -e fastqc_%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

module load fastqc/0.10.1

cd ${and}/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files

fastqc *.fastq.gz --outdir ${and}/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Pcli_fastq_files/
```

Then run multiqc in each directory by just running `multiqc .` (installed it locally on Pegasus scratch space and to PATH variable)

MultiQC reports can be found in the [temperaturevariability2023 GitHub repositorry](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/bioinformatics/AS_pipeline).

## 3) Download Genome [*Acropora cervicornis*](https://usegalaxy.org/u/skitch/h/acervicornis-genome) 

Obtained from [Baums lab](http://baumslab.org/research/data/) with permission from Dr. Sheila Kitchen. Using Version v1.0_171209

Genome file: `Acerv_assembly_v1.0_171209.fasta`

GFF file: `Acerv_assembly_v1.0.gff3`

Protein file: `Acerv_assembly_v1.0.protein.fa`

Transcript file: `Acerv_assembly_v1.0.mRNA.fa`

## 4) Index Acer Genome + Align Experiment Reads to Indexed Genome Using HISAT2 

```{bash}
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
```

## 5) Update .gff3 File for StringTie

I had to manually edit the .gff3 file to make sure the Transcript_ID and Parent_ID were labelled correctly so the IDs match up to the alignment files. R code to edit gff3 file:

```{r}
# Title: A. cervicornis GFF adjustments
# Project: Sedimentation RNA-Seq / Mcap 2020
# Author: J. Ashey / A. Huffmyer -> Allyson DeMerlis
# Date: 06/28/2023

# Need to do some acerv gff adjustments so it can run properly for alignment.

#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="Downloads/Galaxy1-[Acerv_assembly_v1.0.gff3].gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 
Acerv.gff$transcript_id <- gsub("Parent=", "", Acerv.gff$transcript_id) #remove Parent=
head(Acerv.gff)

# Create Parent ID 
Acerv.gff$parent_id <- sub(".*Parent=", "", Acerv.gff$gene)
Acerv.gff$parent_id <- sub(";.*", "", Acerv.gff$parent_id)
Acerv.gff$parent_id <- gsub("ID=", "", Acerv.gff$parent_id) #remove ID= 
head(Acerv.gff)

Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id, ";gene_id=", Acerv.gff$parent_id),  paste0(gene)))
head(Acerv.gff)

Acerv.gff<-Acerv.gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

#save file
write.table(Acerv.gff, file="~/Downloads/Acerv.GFFannotations.fixed_transcript_take3.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```

## 6) Perform Gene Counts with StringTie 

First assemble and estimate reads using aligned .bam files and updated .gff3 file.

```{bash}
#!/bin/bash
#BSUB -J stringtie2_postHISAT2
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie2_postHISAT2%J.out
#BSUB -e stringtie2_postHISAT2%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/AS_pipeline/3_trimmed_fastq_files/Acer_aligned_bam_files"

data=($(ls *.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -p 8 -e -B -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
echo "StringTie assembly for seq file ${i}" $(date) ; \
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

Then, merge StringTie results and obtain gene counts matrix. 

```{bash}
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
```


