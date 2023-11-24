# Bioinformatics pipeline for *P.clivosa*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

Last updated: 2023-11-24

Updates:
1. Changed dds formulas to be ~batch + condition
2. Updated vsd WGCNA to be blind=FALSE
3. change lpv tables so they use FDR-adjusted p-value instead of raw p-value

## Results

Code to create these graphs is from [this R file](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Pcli_Rmd/Pcli_deseq2.R). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment_timepoint)`

Note: in the R script, Treatment_timepoint is actually called 'group' and it is a concatenated factor of Treatment and time_point (i.e., control_Day_29)

countData pre-filtered to remove low-count genes:

```{r}
# how many genes we have total?
nrow(counts) #59947
ncol(counts) #48 samples

# filtering out low-count genes
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData) #34300
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

<img width="751" alt="Screen Shot 2023-08-24 at 1 33 34 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/55cde070-2219-4c40-8c99-5cd81cb33b14">

So because of this, four outliers are removed.

```{r}
outs=c(5,23,39,46) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
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

### 1) PCoA
How many good PCs are there? Look for the number of black points above the line of red crosses (random model). 

<img width="683" alt="Screen Shot 2023-08-24 at 1 34 30 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7aa03c0e-5821-4981-8cf5-a8e1947b5433">

Now plot the PCoA by treatment and time point.

<img width="1118" alt="Screen Shot 2023-09-20 at 12 58 46 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/5d0311c1-4da9-4f54-b872-36dc41a666e2">

Axis 1 is 14.7% of the variance, and axis 2 is 6.0% of the variance.

### 2) PERMANOVA

```{r}
ad=adonis2(t(vsd)~Genotype + time_point*Treatment,data=conditions,method="manhattan",permutations=1e6)
ad
```

<img width="631" alt="Screen Shot 2023-11-24 at 3 21 46 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a0c38276-90b6-41d0-aa2c-8d3531643ec5">

### 3) DESeq2

```{r}
# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)
# treatment
treatment_time=results(dds) 
summary(treatment_time) 
degs_treatment_time=row.names(treatment_time)[treatment_time$padj<0.1 & !(is.na(treatment_time$padj))]
resultsNames(dds)
```

<img width="357" alt="Screen Shot 2023-11-24 at 3 27 50 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/fcc599a6-51eb-4e08-a3fe-89807af861e3">


Specific contrasts:

<img width="344" alt="Screen Shot 2023-08-24 at 1 38 31 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d67dbf56-51fe-44a0-99ba-89ab561ee921">


<img width="349" alt="Screen Shot 2023-08-24 at 1 39 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/85b75332-f670-4ddc-a63a-075e888cd756">


<img width="351" alt="Screen Shot 2023-08-24 at 1 39 34 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/874c65d7-6e61-46b6-8410-b5a95b8747f8">


<img width="341" alt="Screen Shot 2023-08-24 at 1 39 47 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/46c78b3a-5db3-41b0-a31b-a470732df6dd">


### 4) Venn diagram of number of DGEs (based on FDR p-adjusted value < 0.1)

<img width="502" alt="Screen Shot 2023-11-24 at 3 32 48 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/29415867-88a9-41de-9308-bf91992f3dc7">


### 5) PCAs

<img width="622" alt="Screen Shot 2023-08-29 at 4 47 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d0c77eaf-88ff-406d-84c4-1517a15cc785">

<img width="623" alt="Screen Shot 2023-08-29 at 4 47 55 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/84c6fe60-7919-4297-8756-55798456010a">

<img width="618" alt="Screen Shot 2023-09-20 at 1 07 52 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/931961ca-3ac8-4071-a1cf-6732035bc210">

<img width="622" alt="Screen Shot 2023-09-20 at 1 08 09 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/028043e8-dd07-417c-9ee2-8717aab7a4b7">

<img width="627" alt="Screen Shot 2023-09-20 at 1 08 26 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/2182201f-2302-4eec-9124-ef503b4d4c25">

<img width="618" alt="Screen Shot 2023-09-20 at 1 08 45 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/bac6761f-9286-4574-ad9e-f4b60cd9d264">

<img width="621" alt="Screen Shot 2023-09-20 at 1 09 02 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/00150278-6da9-408c-9faa-00b0ffd33268">

<img width="626" alt="Screen Shot 2023-09-20 at 1 09 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d326aa26-6eda-49b2-a10e-5b7c1dec0462">




