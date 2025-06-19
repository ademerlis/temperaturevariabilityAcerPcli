# Bioinformatics pipeline for *Acropora cervicornis*

Last updated: 2025-06-19

## Introduction/Background

I am converting the [Tag-seq](https://github.com/z0on/tag-based_RNAseq) pipeline from Dr. Mikhail Matz at UT Austin from perl(.pl)/python(.py) scripts to bash(.sh) scripts that I can use on the University of Miami High-Performance Computing (HPC) supercomputer, [Pegasus](https://acs-docs.readthedocs.io/pegasus/README.html).

[What is Tag-Seq versus RNA-Seq?](https://pubmed.ncbi.nlm.nih.gov/34674175/#:~:text=3'%2DTag%20RNA%2Dseq,gene%20expression%20by%20tag%20abundance.)

[Meyer et al. 2011 - first publication to use Tag-Seq on a coral species](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2011.05205.x)

Why is it important to follow this exact pipeline for Tag-Seq / 3' RNA-seq / Quant-Seq? Because the 3' method adds a degenerate header at the beginning of each sequence that must be removed. It is also important to remove sequences that share the same degenerate header and the first 20 bases of the sequence. These are duplicate reads and if they are not properly removed, they could conflate the results (by having lots of reads for one gene that are actually just duplicates, not a reflection of gene expression differences). Lastly, it is important to remove the polyA tails and low-quality reads from the sample. [This powerpoint](https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course?preview=%2F103678161%2F225053569%2Ftagseq_for_rnaseqcourse.pdf) has great visuals for this process.

The second point, which is stated twice in the above PowerPoint, is that "pseudo-aligners don't seem to perform as well with Tag-Seq data". Now, I can't find a reference for this, but because this PowerPoint is from the UT Austin Gene Sequencing Facility that my samples were sequenced at, I am going to trust their suggestion. Which means that tools like Salmon are off the table for the quantification step (even though they are touted as being the best and most rapid tool right now).

I am following any code updates made by [Dr. Michael Studivan](https://github.com/mstudiva/tag-based_RNAseq) at NOAA AOML, as he has worked on these codes more recently with our coral species of interest (and has added great notes into his README files).

**The pipeline is as follows:**
0. Download and concatenate reads from Illumina BaseSpace project ([launcher_creator.py](https://github.com/mstudiva/tag-based_RNAseq/blob/master/launcher_creator.py))
1. Run FastQC on raw reads
2. Count number of raw reads ([countreads.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads.pl))
3. Cutadapt to remove adaptors and low-quality reads
4. Download and format reference transcriptomes with bowtie2
5. Concatenate host and symbiont transcriptome
6. Map reads to host/symbiont transcriptomes with bowtie2
7. Count number of mapped reads for mapping efficiency
8. Generate read counts per gene
9. Download gene count files

My adaptations to above pipeline:
1. I skipped the download/concatenate script because I had already downloaded my .fastq.gz files from Basespace. Also I don't have duplicate samples so I didn't need to concatenate any sample files.
2. I previously downloaded [TrimGalore](https://github.com/FelixKrueger/TrimGalore), which is a wrapper for CutAdapt and FastQC. I'm going to try to run the cutadapt script through that since I already have the TrimGalore program locally installed to my Pegasus project space and I have gotten it to work before. The cutadapt version I have installed through TrimGalore(v0.6.10) is v4.4.


## Raw Sequences to Counts Matrix for Importing into R

**Pipeline**: FastQC -> [countreads.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/countreads.pl) -> TrimGalore -> [countreads_trim.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads_trim.pl) -> FastQC -> download and format reference genome or transcriptome -> bowtie2 for index and alignment -> SAMtools for generating counts matrix -> DESeq2

Programs I downloaded locally onto my HPC environment:
- multiQC (version 1.14)
- TrimGalore (version 0.6.10)
- Bowtie2 (version 2.5.2)

Programs I loaded in Pegasus environment that were already installed:
- fastQC (version 0.10.1)
- samtools (version 1.3)

The folder **1_bioinformatics** contains the scripts I used to process the *A. cervicornis* sequencing data on the UM Supercomputer, Pegasus.

### 1. FastQC Raw Reads

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

cd ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/

fastqc *.fastq.gz --outdir ${and}/Ch2_temperaturevariability2023/1_fastq_rawreads/fastqc_files_rawsequences
```

Then you can run multiqc directly in the command line.
```{bash}
cd ${and}/fastqc/
multiqc .
```

The multiqc reports are located in a subfolder on GitHub.

### 2. Count Raw Reads

Made this into a .sh script for Pegasus, but this was originally [countreads.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads.pl). 

```{bash}
#!/bin/bash
#BSUB -J countrawreads
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o countrawreads.out
#BSUB -e countrawreads.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in a bunch of fastq files

#specify variables and paths

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

output_file="countreads_results.txt"

glob=".fastq.gz"
if [ ! -z "$1" ]; then
    glob="$1"
fi

fqs=(*$glob)
for f in "${fqs[@]}"; do
    gunzip -c "$f" > "temp.fastq"  # Decompress the file to a temporary file
    nrd=$(cat "temp.fastq" | wc -l)
    nrd=$((nrd / 4))
    echo "$f	$nrd"
    echo "$f	$nrd" >> "$output_file"  # Append the results to the output file
    rm "temp.fastq"  # Remove the temporary file
done

echo "Results have been saved to $output_file"
```

### 3. Trim Reads

Then I got fancy when I realized I could make a loop within a script that made a bunch of jobs at once, so it could run a bunch of the samples at the same time on Pegasus. 

The parameters I used in cutadapt are specifically from Michael Studivan's code, because it is exactly the degenerate headers and other information that is generated during 3' Tag-Seq and needs to specifically be filtered out. 

This script requires uploading the perl script [tagseq_clipper.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagseq_clipper.pl) to your scripts directory in Pegasus.

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads"

data=($(ls *.fastq.gz))

for samp in "${data[@]}" ; do \

#build script
echo "making cutadapt script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_trim
#BSUB -e ${projdir}/logs/${samp}_trim.err
#BSUB -o ${projdir}/logs/${samp}_trim.out
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/1_fastq_rawreads\"

echo \"Gunzipping file  ${samp} then removing adaptors and low quality reads...\"

gunzip -c ${samp} > ${samp}_temp.fastq
${and}/Ch2_temperaturevariability2023/0_scripts/tagseq_clipper.pl ${samp}_temp.fastq | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${projdir}/${samp/.fastq.gz/}.trim

" > ${projdir}/${samp}_trim.job

#submit script

bsub < ${projdir}/${samp}_trim.job

done

echo "All done!"
```

### 4. Count Trimmed Reads

This .sh script was adapted from [countreads_trim.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads_trim.pl).


```{bash}
#!/bin/bash
#BSUB -J countreads_trim
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o countreads_trim.out
#BSUB -e countreads_trim.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in trimmed fastq files

#specify variables and paths

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4"

output_file="countreads_results.txt"

# Default file pattern
glob="\.trim"

# Check if an argument is provided
if [ "$1" ]; then
    glob="$1"
fi

# Loop through files matching the pattern
for f in *$glob*; do
    # Count the number of lines in the file
    nrd=$(cat "$f" | wc -l)

    # Divide the line count by 4
    nrd=$((nrd / 4))

    # Print the filename and the calculated number
    echo -e "$f\t$nrd" >> "$output_file"
done

echo "Results have been saved to $output_file"
```

### 5. FastQC Trimmed Reads

Again made a loop script to run a bunch of samples at once. 

```{bash}
#BSUB -u and128@miami.edu

#specify variables and paths

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads"


cd ${projdir}/take_4/trimmed_files

data=($(ls *.trim))

for samp in "${data[@]}" ; do \

#build script
echo "making fastqc script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_fastqc
#BSUB -e ${projdir}/logs/${samp}_fastqc.err
#BSUB -o ${projdir}/logs/${samp}_fastqc.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd ${projdir}/take_4/trimmed_files

module load fastqc/0.10.1
fastqc ${samp} --outdir ${projdir}/take_4/trimmed_files/fastqc_output/
echo \"Fastqc script of $samp submitted\"
" > ${projdir}/take_4/trimmed_files/${samp}_fastqc.job

bsub < ${projdir}/take_4/trimmed_files/${samp}_fastqc.job

done
```

### 6. Bowtie2 for Alignment to Reference Genome

First step is to build the reference genome / transcriptome, which comes from *A.cervicornis* [Locatelli et al. 2024](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-11025-3), and  *Symbiodinium* clade A3 from [Shoguchi et al. 2021](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4857-9). 

[Link to download page](https://marinegenomics.oist.jp/symb/viewer/download?project_id=37)
"symatranscriptome_37.fasta.gz"

**Note:** the .bt2 files from the bowtie-build step need to be in the same directory as the trimmed sequences, and cannot be in a subfolder. Otherwise, the alignment code won't be able to identify them. Also, the .bt2 files all need to have the same name before the extensions (i.e. "Locatelli_Shoguchi_concat" in my example), and that is what you write in the alignment script for it to find the files.


```{bash}
#!/usr/bin/env bash
#BSUB -e bowtie2build_LocatelliShoguchiConcat.err
#BSUB -o bowtie2build_LocatelliShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q general
#BSUB -n 8

workdir="/scratch/projects/and_transcriptomics"

bowtie2-build ${workdir}/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa, \
${workdir}/genomes/Symbiodinium/syma_transcriptome_37.fasta \
Locatelli_Shoguchi_concat
```

Second step is alignment.

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer"

data=($(ls *.trim))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_bowtie2align_LocatelliShoguchi
#BSUB -e ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.err
#BSUB -o ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer\"

bowtie2 --local -U ${samp} -x Locatelli_Shoguchi_concat --un ${samp}.unaligned -k 5 -S ${samp}.sam

" > ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

bsub < ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

done
```

### 7. Generate Counts Matrix Using Samtools

This requires uploading the perl script [samcount.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/samcount.pl) to Pegasus before running and putting it in the necessary folder. 

You also need to make a "two-column tab-delimited table transcriptome_seq2gene.tab giving correspondence between entries in the transcriptome fasta file and genes. In de novo transcriptomes, several fasta contigs may correspond to the same gene (e.g., splice isoforms, or assembly uncertainties)." (from [Dr. Matz](https://github.com/z0on/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt))

To generate these files for *A.cervicornis* and *S.fitti*, run this code on the fasta files:

```{bash}
# making seq2iso.tab files
grep ">" Acropora_cervicornis.mrna-transcripts.fa | perl -pe 's/>FUN_(\d+)(\S+)\s.+/FUN_$1$2\t FUN_$1/'>Acervicornis_seq2iso.tab
grep ">" syma_transcriptome_37.fasta | perl -pe 's/>comp(\d+)(\S+)\s.+/comp$1$2\t comp$1/'>Symbiodinium_seq2iso.tab

# create combo file

cat Acer/Locatelli_2023/Acer_Genome/Acervicornis_seq2iso.tab Symbiodinium/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir=

cd "/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch2_temperaturevariability2023/3_genecounts/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch2_temperaturevariability2023/2_trimmed_reads/take_4/trimmed_files/Acer/bowtie2align_LocatelliShoguchi/sam_files\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

bsub < ${and}/Ch2_temperaturevariability2023/3_genecounts/${samp}_samcounts.job

done
```


Then, run these lines of code directly in the command line:

(requires uploading [expressiom_compiler.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/expression_compiler.pl))

```{bash}
perl expression_compiler.pl *.counts > allcounts.txt

# let's remove those annoying chains of extensions from sample names
cat allcounts.txt | perl -pe 's/\_trimmed\.fastq\.gz\.sam\.counts//g'> counts.txt

#and rename FUN -> Acropora and comp -> symbiodinium
sed -i 's/FUN/Acropora/g' counts.txt
sed -i 's/comp/Symbiodinium/g' counts.txt
```

Then, use scp to move the counts.txt file to your local drive.

## Differential Gene Expression Results

Code to create these graphs is from [this R file](https://github.com/ademerlis/temperaturevariabilityAcerPcli/blob/main/Gene%20Expression/Acervicornis/2_DESeq2_host.Rmd). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)`

Treatment has three levels:
1. Initial (day 0)
2. Untreated (control day 29)
3. Treated (variable treated day 29)


### Identifying outliers

Using the R package "biobase", I ran arrayQualityMetrics, which generates a report that visualizes which samples can be considered outliers.

![Screen Shot 2023-12-20 at 3 38 43 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/cbda63be-6894-4d2b-a2b4-abeebdcbbb90)


### PCoA plot

![Screen Shot 2023-12-20 at 3 58 10 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/23cb4d18-e89e-47c8-9c5d-54419c0de6fc)

### PERMANOVA

![Screen Shot 2023-12-20 at 4 01 12 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/dfcba9c5-c521-498b-bea3-22bacb3a8976)

### Significant DGEs

![Screen Shot 2023-12-20 at 4 12 29 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/6319c22e-ac59-41ea-a42e-da18eb06e4ef)

![Screen Shot 2023-12-20 at 4 12 59 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/aa7754a7-eb10-4518-a6d7-ae532304fc1c)

![Screen Shot 2023-12-20 at 4 14 19 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/08c140ab-1d8a-43ae-b658-ae7080981d4b)

![Screen Shot 2023-12-20 at 4 15 16 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/8282e08f-89d5-4ff3-a128-3da2f2f35095)

