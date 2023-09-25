# Bioinformatics pipeline for *A.cervicornis*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

## Results

Code to create these graphs is from [this R file](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Acer_Rmd/Acer_deseq2.R). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment_timepoint)`

Note: in the R script, Treatment_timepoint is actually called 'group' and it is a concatenated factor of Treatment and time_point (i.e., control_Day_29)

countData pre-filtered to remove low-count genes:

```{r}
# how many genes we have total?
nrow(counts) #57358
ncol(counts) #48 samples

# filtering out low-count genes
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData) #47882
ncol(countData) #48
```

First, the package arrayQualityMetrics was used to identify outliers. It generates a report called "index.html", where you can visualize which samples exceed a certain threshold. Then, those are removed from downstream analyses.

```{r}
library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("group"),force=T)
```

Genotype is not included in the intgroup for outliers because it is not the main factor of interest (see [this example](https://github.com/mstudiva/SCTLD-intervention-transcriptomics/blob/main/code/transmission/deseq2_transmission_mcav_host.R)). 

<img width="649" alt="Screen Shot 2023-08-24 at 10 41 55 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/bba2faff-6237-4f13-b4ab-c411f8a11f5f">

So because of this, three outliers are removed.

```{r}
outs=c(20,22,28) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]
```

Then dds model is remade without those outliers. 
```{r}
# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ group + Genotype)
dds$group <- factor(dds$group, levels = c("control_Day_0","control_Day_29","variable_Day_0","variable_Day_29"))
```

### PCoA 

<img width="773" alt="Screen Shot 2023-09-19 at 12 52 02 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/51edd2bf-6ffd-4f65-8f1f-aec5f7db8cb4">

Based on the dds.pcoa table, the % variation explained by each axis is 19.1% for axis 1 and 10.9% for axis 2.

### PERMANOVA

<img width="834" alt="Screen Shot 2023-09-25 at 6 32 32 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/903e9e41-5e3f-4e4a-962a-49e66c500d37">

### DESeq2

```{r}
# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)
# treatment
treatment_time=results(dds) 
summary(treatment_time) 
degs_treatment_time=row.names(treatment_time)[treatment_time$padj<0.1 & !(is.na(treatment_time$padj))]
resultsNames(dds)
```

<img width="368" alt="Screen Shot 2023-09-25 at 6 33 22 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/3bd550b3-c8a0-48e6-977a-5dfe19fb54f2">

<img width="829" alt="Screen Shot 2023-09-25 at 6 33 35 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7f13e4e2-cabb-48e7-b021-995c171c892d">

Specific contrasts:

<img width="604" alt="Screen Shot 2023-08-24 at 10 56 43 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7724d025-2ef0-46d3-aa50-846624cff017">

<img width="567" alt="Screen Shot 2023-08-24 at 10 56 59 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/82877abf-8087-4386-aad5-8d0235f03618">

<img width="608" alt="Screen Shot 2023-08-24 at 10 57 15 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/3907df8b-38f3-4399-9d20-c6833099e85d">

<img width="592" alt="Screen Shot 2023-08-24 at 10 57 30 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a79caefe-0f16-4e58-9b67-a4b1b1039e8e">




