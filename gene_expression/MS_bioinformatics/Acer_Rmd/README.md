# Bioinformatics pipeline for *A.cervicornis*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

## Results

Code to create these graphs is from [this R file]. 

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

### 1. PCoA 

<img width="773" alt="Screen Shot 2023-09-19 at 12 52 02 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/51edd2bf-6ffd-4f65-8f1f-aec5f7db8cb4">

Based on the dds.pcoa table, the % variation explained by each axis is 19.1% for axis 1 and 10.9% for axis 2.

PCoA axes 2 and 3:

![Screen Shot 2023-09-26 at 12 58 33 PM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/cad3af59-4021-4e14-9500-62b97a90bf83)


### 2. PERMANOVA

<img width="834" alt="Screen Shot 2023-09-25 at 6 32 32 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/903e9e41-5e3f-4e4a-962a-49e66c500d37">

### 3. DESeq2

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

### 4. Venn Diagram of number of DGEs (based on FDR p-adjusted value < 0.1)

<img width="641" alt="Screen Shot 2023-09-26 at 10 07 10 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7058529f-0d9a-40ea-a798-24a58901064f">

### 5. PCAs

Plotting everything together, genotype drives clustering.

<img width="570" alt="Screen Shot 2023-09-26 at 10 10 16 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/467c068b-9595-4fbb-882d-5e113b39ce46">

When individual genotypes are plotted:

![Screen Shot 2023-09-26 at 10 35 16 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/12580804-afa5-4802-8593-9d0fc026d811)

![Screen Shot 2023-09-26 at 10 35 33 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/6ebf9435-e180-44e8-8525-0a7bed570e64)

![Screen Shot 2023-09-26 at 10 36 52 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/21260308-f5c2-4f5d-8236-858b67a29599)

Using the Treatment_timepoint grouping variable to plot individual genotypes:

![Screen Shot 2023-09-26 at 10 37 22 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ed0b88db-2e42-4f4c-b9a2-2e339f128f12)

![Screen Shot 2023-09-26 at 10 37 56 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/f1df066d-29d8-43cf-9cdf-4832a7ec39cb)

![Screen Shot 2023-09-26 at 10 38 36 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/13d04de2-5a50-4311-92aa-70f44881fd62)

PC axes 2 and 3 show clearer separation of time point (but not treatment): 

![Screen Shot 2023-09-26 at 10 56 42 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/c5c8f061-161b-420d-9b79-e14db7fa100a)

![Screen Shot 2023-09-26 at 10 56 50 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/4a81b061-c55d-4b21-86d9-d3ac3eff21ab)

![Screen Shot 2023-09-26 at 10 57 22 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/b8c89b5c-f6ae-4542-b0f3-5a56a450d4c4)

### 6. Volcano Plots

![Screen Shot 2023-09-26 at 11 43 42 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a929b4d6-d595-44c5-95a5-2c4feb540719)

![Screen Shot 2023-09-26 at 11 44 23 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/4458fff2-1323-4693-86e4-b00413a84072)

![Screen Shot 2023-09-26 at 11 45 12 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/4c42360e-e8c5-4d3b-840c-17f21d281c31)

![Screen Shot 2023-09-26 at 11 45 46 AM](https://github.com/ademerlis/temperaturevariability2023/assets/56000927/c67bbffc-bdad-4df0-8baf-7c74cedaadd6)

### 7) Common Genes

The main comparison of interest is C0/C29 versus V0/V29. This is because the differential gene expression found in both of these contrasts are solely due to the treatment over time. 

We can find the number of genes that are shared by both these contrasts as well as ones that are unique. The unique ones imply that they are genes specifically changed by the treatment itself (i.e. a DGE in C0/C29 but not in V0/V29 means that the control (untreated) corals resulted in a significant change of expression for that gene but that response was not observed in the variable temperature-treated corals). 

Definition of unique genes: when I look at the lpv table for C0/C29 and V0/V29, they both have the same length = 47,866 isogroups. This is after filtering out all the p-values that are NA. So maybe "unique" isn't the right word, but uniquely significantly differentially expressed is a more specific phrase. 

Venn Diagram of DGEs common in both C0/C29 and V0/V29:


