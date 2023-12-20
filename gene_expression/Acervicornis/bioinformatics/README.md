## Bioinformatics Pipeline (Tag-Seq)

Last updated: 11-09-2023

I am converting the [Tag-seq](https://github.com/z0on/tag-based_RNAseq) pipeline from Dr. Mikhail Matz at UT Austin from perl(.pl)/python(.py) scripts to bash(.sh) scripts that I can use on the University of Miami High-Performance Computing (HPC) supercomputer, [Pegasus](https://acs-docs.readthedocs.io/pegasus/README.html).

[What is Tag-Seq versus RNA-Seq?](https://pubmed.ncbi.nlm.nih.gov/34674175/#:~:text=3'%2DTag%20RNA%2Dseq,gene%20expression%20by%20tag%20abundance.)

[Meyer et al. 2011 - first publication to use Tag-Seq on a coral species](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2011.05205.x)

Why is it important to follow this exact pipeline for Tag-Seq / 3' RNA-seq / Quant-Seq? Because the 3' method adds a degenerate header at the beginning of each sequence that must be removed. It is also important to remove sequences that share the same degenerate header and the first 20 bases of the sequence. These are duplicate reads and if they are not properly removed, they could conflate the results (by having lots of reads for one gene that are actually just duplicates, not a reflection of gene expression differences). Lastly, it is important to remove the polyA tails and low-quality reads from the sample. [This powerpoint](https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course?preview=%2F103678161%2F225053569%2Ftagseq_for_rnaseqcourse.pdf) has great visuals for this process.

The second point, which is stated twice in the above PowerPoint, is that "pseudo-aligners don't seem to perform as well with Tag-Seq data". Now, I can't find a reference for this, but because this PowerPoint is from the UT Austin Gene Sequencing Facility that my samples were sequenced at, I am going to trust their suggestion. Which means that tools like Salmon are off the table for the quantification step (even though they are touted as being the best and most rapid tool right now).

I am following any code updates made by [Dr. Michael Studivan](https://github.com/mstudiva/tag-based_RNAseq) at NOAA AOML, as he has worked on these codes more recently with our coral species of interest (and has added great notes into his README files).

**The pipeline is as follows:**
1. Download and concatenate reads from Illumina BaseSpace project ([launcher_creator.py](https://github.com/mstudiva/tag-based_RNAseq/blob/master/launcher_creator.py))
2. Count number of raw reads ([countreads.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads.pl))
3. Run FastQC on raw reads
4. Cutadapt to remove adaptors and low-quality reads
5. Download and format reference transcriptomes with bowtie2
6. MY ADDITION - based on suggestion from Jenna Dilworth (PhD Candidate from Dr. Carly Kenkel's lab), for the Libro et al. 2013 *Acropora cervicornis* genome, concatenate the *Symbiodinium fitti* genome and run the alignment step with the concatenated file. This is because there is a lot of symbiont read contamination in the supposedly "Acer" genome.
7. Map reads to host/symbiont transcriptomes with bowtie2
8. Count number of mapped reads for mapping efficiency
9. Generate read counts per gene
10. Download gene count files

My adaptations to above pipeline:
1. I skipped the download/concatenate script because I had already downloaded my .fastq.gz files from Basespace. Also I don't have duplicate samples so I didn't need to concatenate any sample files.
2. Mentioned above in #6, but restating here that I am going to concatenate the host and symbiont genomes or transcriptomes for the alignment step.
3. I previoulsy downloaded [TrimGalore](https://github.com/FelixKrueger/TrimGalore), which is a wrapper for CutAdapt and FastQC. I'm going to try to run the cutadapt script through that since I already have the TrimGalore program locally installed to my Pegasus project space and I have gotten it to work before. The cutadapt version I have installed through TrimGalore(v0.6.10) is v4.4 .
