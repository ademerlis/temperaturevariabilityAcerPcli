# Bioinformatics pipeline for *P.clivosa*

To do: add scripts from [Michael's Tag-based_RNAseq pipeline](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt)  that he ran in his HPC.

Last updated: 2023-11-24

Updates:
1. Changed dds formulas to be ~batch + condition
2. Updated vsd WGCNA to be blind=FALSE
3. change lpv tables so they use FDR-adjusted p-value instead of raw p-value
4. Combine control and variable day 0 into one group, called "Initial" (because no treatment has happened yet, so separating them doesn't make sense)


## Results

`dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)`

Treatment has three levels:
1. Initial (day 0)
2. Untreated (control day 29)
3. Treated (variable treated day 29)

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
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)
dds$Treatment <- factor(dds$Treatment, levels = c("Initial", "Untreated", "Treated"))
```

### 1) PCoA

<img width="531" alt="Screen Shot 2023-11-24 at 7 59 04 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/56012fac-03ec-4af5-9a52-e6a53ed5f72b">

First axis explains 14.8% of the variance, axis 2 explains 6% of the variance.

### 2) PERMANOVA

```{r}
ad=adonis2(t(vsd)~Genotype + Treatment,data=conditions,method="manhattan",permutations=1e6)
ad
```

<img width="731" alt="Screen Shot 2023-11-24 at 8 01 54 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/04864f5b-e4dd-46e5-940b-45c3ef13c2d1">


### 3) DESeq2

```{r}
# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)
# treatment
Treatment=results(dds) 
summary(Treatment) 
degs_Treatment=row.names(Treatment)[Treatment$padj<0.1 & !(is.na(Treatment$padj))]
resultsNames(dds)
```

**some notes from DESeq2 when it was running:**
-- note: fitType='parametric', but the dispersion trend was not well captured by the
   function: y = a/x + b, and a local regression fit was automatically substituted.
   specify fitType='local' or 'mean' to avoid this message next time.
final dispersion estimates, fitting model and testing: 6 workers
-- replacing outliers and refitting for 24 genes
-- DESeq argument 'minReplicatesForReplace' = 7 

<img width="349" alt="Screen Shot 2023-11-24 at 8 05 32 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/08c2aad3-a438-4692-9a23-fdb971288e1e">


Specific contrasts:

**How to interpret**: the contrast reported first is the numerator, i.e. Treated (numerator) vs. Untreated (denomenator), the number of DEGs upregulated are differentially upregulated (greater positive LogFoldChange) in the treated group than the untreated group. 

<img width="352" alt="Screen Shot 2023-11-24 at 8 07 39 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/48b3409a-59ea-4dc1-9637-02d38867ee5a">

<img width="351" alt="Screen Shot 2023-11-24 at 8 08 03 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a5213278-4016-4c5b-8dd7-39c8c2f7d335">

<img width="349" alt="Screen Shot 2023-11-24 at 8 10 44 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/134a35b7-ea1c-4577-9fcc-a0c0fa413a19">


### 4) Venn diagram of number of DGEs (based on FDR p-adjusted value < 0.1)

<img width="447" alt="Screen Shot 2023-11-24 at 8 45 22 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ab06e6cf-739a-443c-b74d-232552a7c4d3">

This one was made with VennDiagram and it's ugly and I can't move the labels or make the circles smaller.

This one was made with ggvenn and is marginally better:

<img width="612" alt="Screen Shot 2023-11-24 at 8 46 09 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/e9f995b2-0f1d-4ff7-bf45-45b6564a56fc">


### 5) PCAs

Plotting everything together, genotype drives clustering.

<img width="549" alt="Screen Shot 2023-11-24 at 3 54 07 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/71e4605c-acf4-4037-a68b-537c904c17c9">

When individual genotypes are plotted: 

<img width="625" alt="Screen Shot 2023-11-24 at 9 04 39 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d4daa57c-826e-4725-90fe-9a0551476942">

But with the stat_ellipse function it looks better:

<img width="620" alt="Screen Shot 2023-11-24 at 9 05 10 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/278e48ed-07cd-4720-80de-d6c141a1efbe">

<img width="625" alt="Screen Shot 2023-11-24 at 9 05 36 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/4666ecf5-d8b5-49f5-b2ae-e98834cba0d1">

<img width="623" alt="Screen Shot 2023-11-24 at 9 05 58 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/c5544b57-97ae-45ec-84df-ae8a1f9366fc">

<img width="622" alt="Screen Shot 2023-11-24 at 9 06 27 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/2efc4e87-5aee-4bb4-b715-b94fd8f0d593">

<img width="621" alt="Screen Shot 2023-11-24 at 9 06 48 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/337fd5ba-b500-4b86-b617-a69c2e35da3d">

the PC axes 2 and 3 don't show any patterns (other than genotype again).

<img width="623" alt="Screen Shot 2023-11-24 at 9 08 41 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/446179d9-ef82-42e1-9a4f-a39515cafa54">

### 6) Volcano plots

<img width="585" alt="Screen Shot 2023-11-24 at 9 11 12 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/b5e11831-6a7d-422c-868e-e1587e385151">

<img width="587" alt="Screen Shot 2023-11-24 at 9 12 28 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/e9c29af2-7de5-4648-98e9-257072675998">

<img width="585" alt="Screen Shot 2023-11-24 at 9 13 40 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/17555bb8-a018-4292-a7b4-e1603be8582a">



