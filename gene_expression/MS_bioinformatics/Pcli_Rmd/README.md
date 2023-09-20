# Bioinformatics pipeline for *P.clivosa*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

## Results

Code to create these graphs is from [this R file](https://github.com/ademerlis/temperaturevariability2023/blob/main/gene_expression/MS_bioinformatics/Pcli_Rmd/Pcli_deseq2.R). 

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ group + Genotype)`

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

### 1) Heatmap
To see similarity of samples

<img width="681" alt="Screen Shot 2023-08-24 at 1 34 14 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/5cb8c257-a7fc-442c-812f-09979d91287e">


### 2) PCoA
How many good PCs are there? Look for the number of black points above the line of red crosses (random model). 

<img width="683" alt="Screen Shot 2023-08-24 at 1 34 30 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7aa03c0e-5821-4981-8cf5-a8e1947b5433">


Now plot the PCoA by treatment and time point.

<img width="700" alt="Screen Shot 2023-08-24 at 1 34 43 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/1f327ca2-8f62-475d-9e36-914b9b8c8a22">


Neighbor-joining tree of samples (based on significant PCoA's).

<img width="631" alt="Screen Shot 2023-08-24 at 1 34 54 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/2e85ccc7-da1a-448a-9139-e8da0d251184">


### 3) PERMANOVA for variance in distance matrices

```{r}
ad=adonis2(t(vsd)~time_point*Treatment + Genotype,data=conditions,method="manhattan",permutations=1e6)
ad
```

<img width="591" alt="Screen Shot 2023-08-24 at 1 35 23 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7ebd17c0-6e3c-4df0-9c46-e2c77bc283b8">


Pie chart to show proportion of R2 values per factor driving variance

<img width="610" alt="Screen Shot 2023-08-24 at 1 36 22 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/834c22b9-5c10-431f-a549-1b0df39d7094">


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

<img width="353" alt="Screen Shot 2023-08-24 at 1 37 51 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ea6804b6-e3f8-48bf-8e7b-07a09cc0a158">


Specific contrasts:

<img width="344" alt="Screen Shot 2023-08-24 at 1 38 31 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d67dbf56-51fe-44a0-99ba-89ab561ee921">


<img width="349" alt="Screen Shot 2023-08-24 at 1 39 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/85b75332-f670-4ddc-a63a-075e888cd756">


<img width="351" alt="Screen Shot 2023-08-24 at 1 39 34 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/874c65d7-6e61-46b6-8410-b5a95b8747f8">


<img width="341" alt="Screen Shot 2023-08-24 at 1 39 47 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/46c78b3a-5db3-41b0-a31b-a470732df6dd">



### 5) Density plot for DEGs

<img width="668" alt="Screen Shot 2023-08-24 at 1 42 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/eecfbb5f-9dd9-4949-b0ad-4ee24908568b">


### 6) Venn diagram for DEGs

<img width="605" alt="Screen Shot 2023-08-24 at 1 49 24 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/f4d287ce-ed29-4bea-a509-79722870c896">

### 7) PCAs

<img width="622" alt="Screen Shot 2023-08-29 at 4 47 17 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d0c77eaf-88ff-406d-84c4-1517a15cc785">

<img width="623" alt="Screen Shot 2023-08-29 at 4 47 55 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/84c6fe60-7919-4297-8756-55798456010a">


