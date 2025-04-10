#[![DOI](https://zenodo.org/.svg)](https://zenodo.org/doi)

[Authorea link](https://www.authorea.com/users/906749/articles/1280948-species-level-differences-in-molecular-responses-to-a-thermally-variable-stress-hardening-treatment-for-caribbean-corals)

This repository contains data and analysis scripts for the manuscript:

## Species-level differences in molecular responses to a thermally variable stress-hardening treatment for Caribbean corals
#### **Authors:** Allyson DeMerlis, Michael S. Studivan, Kevin Wong, Nash Soderberg, David Ehrens, Lys M. Isma, Katrina Rosing, Katrina Sophia Cocson, Rowan Thomas, Danielle Dvorkin, Patrick M. Kiel, Joseph D. Unsworth, Martine D’Alessandro, Ana M. Palacio-Castro, Diego Lirman, Andrew C. Baker, Erinn M. Muller, Nikki Traylor-Knowles, Ian C. Enochs
#### **Journal:** _Submitted to Molecular Ecology_ [doi:XXX](http://dx.doi.org/XXX)  

-----

### Description:
These repository contains all data and code used to study the impact of thermal variability as an ex-situ stress-hardening treatment on the physiology and gene expression of _Acropora cervicornis_ and _Pseudodiploria clivosa_.

### Contents:

#### Tank Parameters:
* **tank_temp_figs.Rmd:** has all the code needed for Figure 1. It sources the data from the **Data from LabVIEW** subfolder and saves statistical summaries and plots in the **plots** subfolder.

#### Physiology:
* The complete data file that contains all metadata and buoyant weight mesaurements is the file **metadata.csv**.
* Physiology data analysis is broken up into three subfolders: **Calcification**, **R_intensity**, and **Photosynthetic Efficiency**.
* In the subfolder **Calcification**, **calcification.Rmd** has all the code needed for Supplementary Figure 1 and Table S3 statistics.
* In the subfolder **R_intensity**, **colorscores.Rmd** has all the code needed for Supplementary Figure 2 and Table S4 statistics.
* The subfolder **Photosynthetic Efficiency** has several R files that are important for analysis.
  - **1_importtidyPAMdata.Rmd** includes the function to import raw IPAM data into R, and creates a format file that matches coral fragment metadata to the area of interest (AOI) and YII (photosynthetic efficiency) values to the IPAM image metadata files. The raw files for this can be found in the **ipam_data** subfolder.
  - **2_treatment_stats.Rmd** has code needed for Figure 2 and Table S4 statistics.
  - **3_rapidheatstressassay_fvfm.Rmd** has code needed for Figure 2 and Table S5 statistics.
 
#### Gene Expression:
* The raw sequence .fastq files will be made publicly available on the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) under BioProject PRJNA1196005, upon publication.
* The file **readcounts_rawtrimmedaligned_summary.csv** has a description of the raw, trimmed, and aligned reads for *A. cervicornis* and *P. clivosa*.
* The file **treatment_metadata.csv** has the sample metadata, including date sampled, treatment group, and genotype.
* For each species, there is a corresponding subfolder.
* In **Acervicornis**:
  - **1_bioinformatics** subfolder contains the scripts needed to process raw 3' RNA-Seq sequences on the UM HPC, [Pegasus](https://acs-docs.readthedocs.io/pegasus/README.html), which uses an LSF resource manager. The bioinformatics pipeline outlined in this folder includes: FastQC > Cutadapt > bowtie2 > samtools
  - **2_DESeq2_host.Rmd** R markdown file for analyzing *A. cervicornis* host differential gene expression as a result of the variable temperature treatment. This has the code needed for Figures 3 and 4, as well as Supplementary Tables 8, 9, and 11.
  - **2_DESeq2_symbiont.Rmd** R markdown file for analyzing *A. cervicornis* symbiont differential gene expression as a result of the variable temperature treatment. This has the code needed for Figures 3 and 4, as well as Supplementary Tables 8, 9, and 12.
  - **2_outlier_detection** subfolder for the *ArrayQualityMetrics* outlier detection method used for removing samples based on gene expression patterns.
  - **3_GO-MWU** subfolder contains the R markdown files and functions necessary for running the specific method of Gene Ontology (GO) enrichment analysis used in this study. This has the code needed for Figures 5 and 6, as well as Supplementary Tables 10, 15, 16, and 17.
  - **results_csv** subfolder contains all the results files for the *A. cervicornis* host and symbiont differential gene expression.
* In **Pclivosa**:

</br>

#### Notes
* This README.md was formatted following [Dr. Ana Palacio's GitHub repository](https://github.com/anampc/Acer_NH4_disease/tree/master).

