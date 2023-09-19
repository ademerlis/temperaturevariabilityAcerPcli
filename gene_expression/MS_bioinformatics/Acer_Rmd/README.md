# Bioinformatics pipelines for *A.cervicornis*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

## Results

1. [Heatmap](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#1-heatmap)
2. [PCoA](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#2-pcoa)
3. [PERMANOVA](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#3-permanova-for-variance-in-distance-matrices)
4. [DESeq2](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#4-deseq2)
5. [Density plot](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#5-density-plot-for-degs)
6. [Venn Diagram of DGEs](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#6-venn-diagram-for-degs)
7. [PCAs](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#7-pca)
8. [KOG-MWU heatmaps](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Acer_Rmd/README.md#8-kog-mwu-heatmaps)
9. [KOG-MWU correlation plots](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Acer_Rmd/README.md#9-kog-mwu-correlation-plots)
10. [Volcano Plots](https://github.com/ademerlis/temperaturevariability2023/tree/main/gene_expression/MS_bioinformatics/Acer_Rmd#10-volcano-plots)
11. WGCNA
12. GO-MWU Analysis of WGCNA Modules 
13. Orthofinder

Code to create these graphs is from [this R file](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Acer_Rmd/Acer_deseq2.R). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ group + Genotype)`

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

the "group" variable for the design contains both "Treatment" (control vs. variable) and "time_point" (Day_0 vs. Day_29).

```{r}
design$Genotype <- as.factor(design$Genotype)
design$Treatment <- as.factor(design$Treatment)
design %>% 
  mutate(time_point = case_when(Experiment.phase == "Pre-treatment" ~ "Day_0",
                                Experiment.phase == "last day of treatment" ~ "Day_29")) -> design
column_to_rownames(design, var="Sample_ID") -> design
design$group <- factor(paste0(design$Treatment, "_", design$time_point))
# reorders fate factor according to "control" vs "treatment" levels
dds$group <- factor(dds$group, levels = c("control_Day_0","control_Day_29","variable_Day_0","variable_Day_29"))
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

### 1) Heatmap
To see similarity of samples

<img width="692" alt="Screen Shot 2023-08-24 at 10 45 03 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/11a824df-3d22-4704-a79a-95affcb2b657">

### 2) PCoA
How many good PCs are there? Look for the number of black points above the line of red crosses (random model). 

<img width="667" alt="Screen Shot 2023-08-24 at 10 46 52 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/184913c8-3b42-42d3-a090-f6d857f8138c">

Now plot the PCoA by treatment and time point.

<img width="773" alt="Screen Shot 2023-09-19 at 12 52 02 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/51edd2bf-6ffd-4f65-8f1f-aec5f7db8cb4">

Based on the dds.pcoa table, the % variation explained by each axis is 19.1% for axis 1 and 10.9% for axis 2.

Neighbor-joining tree of samples (based on significant PCoA's).

<img width="636" alt="Screen Shot 2023-08-24 at 10 48 34 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a9d877ab-7bc1-4515-8ef5-2b60d0dca33f">

### 3) PERMANOVA for variance in distance matrices

```{r}
ad=adonis2(t(vsd)~time_point*Treatment + Genotype,data=conditions,method="manhattan",permutations=1e6)
ad
summary(ad)
```
<img width="735" alt="Screen Shot 2023-08-24 at 10 51 22 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/c665756c-7148-4cf0-9d4b-22d8bae042f9">

Pie chart to show proportion of R2 values per factor driving variance

<img width="423" alt="Screen Shot 2023-08-24 at 10 52 07 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/637c63d6-590c-42f2-8b70-b66cacf82618">

### 4) DESeq2

```{r}
# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)
# treatment
treatment_time=results(dds) 
summary(treatment_time) 
degs_treatment_time=row.names(treatment_time)[treatment_time$padj<0.1 & !(is.na(treatment_time$padj))]
resultsNames(dds)
```

<img width="348" alt="Screen Shot 2023-08-24 at 10 55 33 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/e90ca043-d909-4b34-905b-2d7af2d8e141">

<img width="565" alt="Screen Shot 2023-08-24 at 10 56 15 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/41649c31-9fec-4aea-83d6-804a22bb7a7e">

Specific contrasts:

<img width="604" alt="Screen Shot 2023-08-24 at 10 56 43 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7724d025-2ef0-46d3-aa50-846624cff017">

<img width="567" alt="Screen Shot 2023-08-24 at 10 56 59 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/82877abf-8087-4386-aad5-8d0235f03618">

<img width="608" alt="Screen Shot 2023-08-24 at 10 57 15 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/3907df8b-38f3-4399-9d20-c6833099e85d">

<img width="592" alt="Screen Shot 2023-08-24 at 10 57 30 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a79caefe-0f16-4e58-9b67-a4b1b1039e8e">

### 5) Density plot for DEGs

<img width="462" alt="Screen Shot 2023-08-24 at 10 58 13 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/664f279e-f805-479c-8998-62c8774bd7f8">

### 6) Venn diagram for DEGs

<img width="573" alt="Screen Shot 2023-08-24 at 10 58 58 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/50f945b9-0b34-4601-afe3-c7adac07bb01">

### 7) PCA

<img width="625" alt="Screen Shot 2023-08-24 at 2 39 15 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7d7729dc-09a2-4271-a42d-f3576d63de99">

<img width="628" alt="Screen Shot 2023-08-24 at 2 39 44 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/9f8b6e55-a573-4252-9c6d-84ddddc1710f">

Separating each genotype out:

1. SI-C

<img width="629" alt="Screen Shot 2023-09-16 at 7 10 27 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/e72fa911-130e-44cd-a478-9af22be82e36">

<img width="629" alt="Screen Shot 2023-09-16 at 7 10 57 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/1abd1c9f-d34f-4bee-a6bb-b180f2f76f31">

2. MB-B

<img width="623" alt="Screen Shot 2023-09-16 at 7 12 15 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/16286689-f9e3-4790-a9cb-b11223c996c5">

<img width="629" alt="Screen Shot 2023-09-16 at 7 12 38 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/67760487-2401-4e98-8ec8-b8139426ae52">

3. BC-8b

<img width="623" alt="Screen Shot 2023-09-16 at 7 13 48 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/83b12325-7d5d-4621-b724-3d1f1d290e6c">

<img width="624" alt="Screen Shot 2023-09-16 at 7 14 25 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/908c3d1f-e72e-45c0-b201-4bd5ccb1940f">

Too few samples on Day 29 to create an ellipse.

Then I looked at PC axes 2 and 3, and there was still a separation of genotype, but also a clear separation of time point.

<img width="622" alt="Screen Shot 2023-09-16 at 7 51 19 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/914f1a47-729f-4316-8ea7-bbc904893d86">

<img width="631" alt="Screen Shot 2023-09-16 at 7 51 35 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/5b8b671c-4b0c-4cd1-bbcb-24610ecbcd46">

### 8) KOG-MWU Heatmaps

1. Enrichment of KOG (euKaryotic Orthologous Groups) terms of DGEs with correlations based on **log-transformed p-values**. Note that the scale of the heatmap is not p-values, but delta-ranks, which are generated based on Mann-Whitney U tests (where the MWU comes from).
    
<img width="605" alt="Screen Shot 2023-09-16 at 8 03 29 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/0d5300c4-4f92-4d73-9543-f7bcb2f7cc9e">

2. Enrichment of KOG (euKaryotic Orthologous Groups) terms of DGEs with correlations based on **fold change**. Note that the scale of the heatmap is not fold change, but delta-ranks, which are generated based on Mann-Whitney U tests (where the MWU comes from).

<img width="579" alt="Screen Shot 2023-09-16 at 8 06 29 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/071364f7-a0b5-465f-b804-b74b7b4f64f1">


### 9) KOG-MWU Correlation plots

1. Correlation plots for pairwise comparisons of treatment_time based on KOG enrichment for **log-transformed p-values**.
   
<img width="692" alt="Screen Shot 2023-09-16 at 8 08 02 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/bd2c6b13-afa6-4b0e-a1c0-8462e09e312f">

2. Correlation plots for pairwise comparisons of treatment_time based on KOG enrichment for **fold change**.
   
<img width="689" alt="Screen Shot 2023-09-16 at 8 08 52 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/0f1604ff-079f-466d-9766-7e9ddaf46b10">

### 10) Volcano Plots

NOTE: for all volcano plots, the y-axis is p-adjusted value, I just didn't change the label to say that and got lazy.

<img width="1227" alt="Screen Shot 2023-09-17 at 6 25 24 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/e46aba63-8ebd-4811-b40e-87efae869026">

<img width="1235" alt="Screen Shot 2023-09-17 at 6 29 05 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/19419923-7fe5-47c2-9ff3-a09ba3909f0f">

<img width="1227" alt="Screen Shot 2023-09-17 at 6 31 30 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/0a04e1a4-d2c5-4b4c-9a59-457ffd882d5e">

<img width="1230" alt="Screen Shot 2023-09-17 at 6 32 14 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/aa450d8c-31c5-46c1-a0d9-d606d3d68a9f">

### 11) WGCNA

Code for WGCNA is found [here](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Acer_Rmd/WGCNA/wgcna_Acer.R). 

First part of WGCNA is outlier detection using goodSamplesGenes. No outliers detected here. Then plot the log-transformed variance-stabilized gene expression across all samples and look for any obvious deviations (there are none). 

<img width="736" alt="Screen Shot 2023-09-19 at 1 03 49 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/86d46996-0de1-45e8-96b8-c525a481d8cc">

Then plot a sample dendrogram to look for outliers.

<img width="868" alt="Screen Shot 2023-09-19 at 1 05 58 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d324cbe2-2e92-408b-aa83-03568f03d115">

One outlier was found, which was then removed using this code:

```{r}
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datt=datt[!remove.samples,]
traits=traits[!remove.samples,] #1 sample removed
```

Next is the soft threshold. 

<img width="769" alt="Screen Shot 2023-09-19 at 1 06 55 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/2f1f7521-cf54-483a-9f9c-a9bf531fe938">

The goal is to get something that corresponds to an R^2 cutoff of 0.90, but none of mine reached that. The closest is a sft power of 21 which had an sft R^2 of 0.897. 

Then run adjacency -> TOM -> dissTOM to make a gene tree. Then make modules (started with minimum module size of 30), and then calculate eigengenes from there. 

Do 2 passes of module generation, cluster dendrogram, and module correlations. In the first pass, don't set a threshold for module eigengene dissimilarity because we want to see how the module-trait heatmap is generated, then decide which modules look really similar and can be merged. 

<img width="769" alt="Screen Shot 2023-09-19 at 1 11 36 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a4cb57b8-a840-4f61-8620-65f261a40f56">

In the second pass, the ME dissimilarity threshold was set to 0.4 because that's what Michael used and it also looked like a good cut-off point. 

<img width="753" alt="Screen Shot 2023-09-19 at 1 11 55 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/b940d281-d63b-4f71-a834-0e878aa33ad9">

Here are the unmerged vs merged modules in the dendrogram.

<img width="759" alt="Screen Shot 2023-09-19 at 1 12 14 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/efbde41c-6be0-403d-97f3-af624a826138">

Next, we correlate the modules to all the "traits" (conditions in the DESeq2 formula). 

<img width="728" alt="Screen Shot 2023-09-19 at 1 14 23 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/fc26f843-c583-483d-8937-b0cafe03049a">

Module size barplot: 

<img width="656" alt="Screen Shot 2023-09-18 at 3 20 00 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/93ce5182-6d94-495c-8488-4609561d10e7">







