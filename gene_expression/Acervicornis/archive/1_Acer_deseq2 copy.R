#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("DESeq2",dependencies=T)
#BiocManager::install("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# BiocManager::install("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")


#### DATA IMPORT ####
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)
library(tidyverse)

#read in counts
counts = read.delim("final_counts_matrix.txt")

#select only Acropora genes, remove Symbiodinium
counts <- counts %>% filter(!grepl("^Symbiodinium[0-9]+$", X))

column_to_rownames(counts, var ="X") -> counts
counts %>% 
  select(!X.1) -> counts

#two of the sample IDs are in all caps, need to change it to be "Acer"
counts <- counts %>% dplyr::rename(Acer.112 = ACER.112)
counts <- counts %>% dplyr::rename(Acer.124 = ACER.124)
 
# how many genes we have total?
nrow(counts) #30122
ncol(counts) #48 samples

# how does the data look? 
head(counts)

colnames(counts) # sample IDs (i.e. Acer.005). You need to keep track of column names and row names of data frames throughout this analysis because they MUST MATCH 
rownames(counts) #gene names (but not actual gene names, arbitrary gene names like Acropora000001. To get "actual" gene names, need annotation files)

# filtering out low-count genes
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData) #25003
ncol(countData) #48
#write.csv(countData, file = "Acer_countdata.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna) #15526
ncol(counts4wgcna) #48
#write.csv(counts4wgcna, file="Acer_counts4wgcna.csv")

colnames(counts4wgcna) #still same sample ID names "Acer.005" etc

# importing a design .csv file
design = read.csv("../treatment_metadata.csv", head=TRUE)

design %>% 
  filter(Species == "Acer") -> design

design %>% 
  select(!Species) -> design

column_to_rownames(design, var="Sample_ID") -> design
design$Genotype <- as.factor(design$Genotype)
design$Genotype <- factor(gsub("-", "_", design$Genotype)) #DESeq2 does not like hyphens in factor names
design$Treatment <- as.factor(design$Treatment)

#### FULL MODEL DESIGN (Genotype + Treatment) and OUTLIERS ####

#when making dds formula, it is CRITICAL that you put the right order of variables in the design. The design indicates how to model the samples, 
#(here: design = ~batch + condition), 
#that we want to measure the effect of the condition, controlling for batch differences. The two factor variables batch and condition should be columns of coldata.

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

colnames(Vsd)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("Treatment"),force=T) #Genotype is not included as an intgroup because it is not the main factors of interest.
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested.
#I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. 
#Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(20,22,28,36) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="initial_fullddsdesigncountsVsdcountsWGCNA.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
colnames(vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,1], design[,4],sep="_")
snames #i.e. 

# renames the column names
colnames(vsd)=snames

save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ Genotype + Treatment)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=FALSE) #blind=TRUE is the default, and it is a fully unsupervised transformation. However, the creator of DESeq2,
#Michael Love, recommends using blind=FALSE for downstream analyses because when transforming data, the full use of the design information should be made. If many genes have
#large differences in counts due to experimental design, then blind=FALSE will account for that.

head(vsd.wg)
colnames(vsd.wg)=snames
colnames(vsd.wg)
colnames(vsd)
save(vsd.wg,design,file="data4wgcna.RData")

#### DESEQ ####

# with multi-factor, multi-level design
load("Rdata_files/initial_fullddsdesigncountsVsdcountsWGCNA.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
#dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)
rownames(design)
colnames(countData)
dds=DESeq(dds, parallel=TRUE)

# saving all models
save(dds,file="Rdata_files/realModels_Acer.RData")

#### DEGs and CONTRASTS ####

load("Rdata_files/realModels_Acer.RData")
library(DESeq2)
library(tidyverse)

read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                          sep = "\t",
                          quote="", fill=FALSE) %>%
  rename(gene = V1,
         annot = V2) -> iso2geneName

# treatment
Treatment=results(dds) 
summary(Treatment, alpha = 0.05) 
resultsNames(dds)
#primary ones of interest: "Treatment_Treated_vs_Initial",  "Treatment_Untreated_vs_Initial",  "Treatment_Treated_vs_Untreated"

degs_treatment=row.names(Treatment)[Treatment$padj<0.05 & !(is.na(Treatment$padj))]
length(degs_treatment) #7272 genes

# Treated vs. Untreated
Treatment_Treated_vs_Untreated=results(dds,contrast=c("Treatment","Treated","Untreated"))
summary(Treatment_Treated_vs_Untreated, alpha = 0.05)
degs_Treatment_Treated_vs_Untreated=row.names(Treatment_Treated_vs_Untreated)[Treatment_Treated_vs_Untreated$padj<0.05 & !(is.na(Treatment_Treated_vs_Untreated$padj))]
length(degs_Treatment_Treated_vs_Untreated) #1587

Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange > 1) %>% 
  nrow() #200
Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange < -1) %>% 
  nrow() #207

Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("TreatedvUntreated_annotDGEs.csv")

Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  mutate(Contrast = "Treated vs Untreated") -> treated_vs_untreated

# Treated vs. Initial
Treatment_Treated_vs_Initial=results(dds,contrast=c("Treatment","Treated","Initial"))
summary(Treatment_Treated_vs_Initial, alpha = 0.05)
degs_Treatment_Treated_vs_Initial=row.names(Treatment_Treated_vs_Initial)[Treatment_Treated_vs_Initial$padj<0.05 & !(is.na(Treatment_Treated_vs_Initial$padj))]
length(degs_Treatment_Treated_vs_Initial) #5733

Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange > 1) %>% 
  nrow() #489
Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange < -1) %>% 
  nrow() #565

Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("TreatedvInitial_annotDGEs.csv")

Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  mutate(Contrast = "Treated vs Initial") -> treated_vs_initial


# Untreated vs. Initial
Treatment_Untreated_vs_Initial=results(dds,contrast=c("Treatment","Untreated","Initial"))
summary(Treatment_Untreated_vs_Initial, alpha = 0.05)
degs_Treatment_Untreated_vs_Initial=row.names(Treatment_Untreated_vs_Initial)[Treatment_Untreated_vs_Initial$padj<0.05 & !(is.na(Treatment_Untreated_vs_Initial$padj))]
length(degs_Treatment_Untreated_vs_Initial) #7272


Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange > 1) %>% 
  nrow() #1050
Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj<0.05 & log2FoldChange < -1) %>% 
  nrow() #943

Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>%
  write_csv("UntreatedvInitial_annotDGEs.csv")

Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  mutate(Contrast = "Untreated vs Initial") -> untreated_vs_initial

rbind(treated_vs_untreated, treated_vs_initial, untreated_vs_initial) %>% 
  write_csv("all_contrasts_annotDEGs.csv")


save(Treatment,Treatment_Untreated_vs_Initial,Treatment_Treated_vs_Initial,Treatment_Treated_vs_Untreated, degs_Treatment_Untreated_vs_Initial,degs_Treatment_Treated_vs_Initial, degs_Treatment_Treated_vs_Untreated, file="Rdata_files/pvals.RData")

## lfcThreshold=1
# Treated vs. Untreated
Treatment_Treated_vs_Untreated_LFC1=results(dds,contrast=c("Treatment","Treated","Untreated"), lfcThreshold=1)
summary(Treatment_Treated_vs_Untreated_LFC1, alpha = 0.05)
degs_Treatment_Treated_vs_Untreated_LFC1=row.names(Treatment_Treated_vs_Untreated_LFC1)[Treatment_Treated_vs_Untreated_LFC1$padj<0.05 & !(is.na(Treatment_Treated_vs_Untreated_LFC1$padj))]
length(degs_Treatment_Treated_vs_Untreated_LFC1) #10

# Treated vs. Initial
Treatment_Treated_vs_Initial_LFC1=results(dds,contrast=c("Treatment","Treated","Initial"), lfcThreshold=1)
summary(Treatment_Treated_vs_Initial_LFC1, alpha = 0.05)
degs_Treatment_Treated_vs_Initial_LFC1=row.names(Treatment_Treated_vs_Initial_LFC1)[Treatment_Treated_vs_Initial_LFC1$padj<0.05 & !(is.na(Treatment_Treated_vs_Initial_LFC1$padj))]
length(degs_Treatment_Treated_vs_Initial_LFC1) #35

# Untreated vs. Initial
Treatment_Untreated_vs_Initial_LFC1=results(dds,contrast=c("Treatment","Untreated","Initial"), lfcThreshold=1)
summary(Treatment_Untreated_vs_Initial_LFC1, alpha = 0.05)
degs_Treatment_Untreated_vs_Initial_LFC1=row.names(Treatment_Untreated_vs_Initial_LFC1)[Treatment_Untreated_vs_Initial_LFC1$padj<0.05 & !(is.na(Treatment_Untreated_vs_Initial_LFC1$padj))]
length(degs_Treatment_Untreated_vs_Initial_LFC1) #95


#### VENN DIAGRAMS ####

load("RData_files/pvals.RData")
library(DESeq2)

pairwise=list("Untreated vs. Initial"=degs_Treatment_Untreated_vs_Initial, "Treated vs. Initial"=degs_Treatment_Treated_vs_Initial,"Treated vs. Untreated"=degs_Treatment_Treated_vs_Untreated)

library(ggvenn)

ggvenn(pairwise) + 
  scale_fill_manual(values = c("#ca0020", "#0571b0", "#f4a582"))

pairwise_LFC1=list("Untreated vs. Initial"=degs_Treatment_Untreated_vs_Initial_LFC1, "Treated vs. Initial"=degs_Treatment_Treated_vs_Initial_LFC1,"Treated vs. Untreated"=degs_Treatment_Treated_vs_Untreated_LFC1)

ggvenn(pairwise_LFC1) + 
  scale_fill_manual(values = c("#ca0020", "#0571b0", "#f4a582"))

#### GO/KOG EXPORT ####

load("RData_files/realModels_Acer.RData")
load("RData_files/pvals.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# Untreated vs Initial
# log2 fold changes:
source=Treatment_Untreated_vs_Initial[!is.na(Treatment_Untreated_vs_Initial$padj),]
Untreated_vs_Initial.fc=data.frame("gene"=row.names(source))
Untreated_vs_Initial.fc$lfc=source[,"log2FoldChange"]
head(Untreated_vs_Initial.fc)
write.csv(Untreated_vs_Initial.fc,file="Untreated_vs_Initial_fc.csv",row.names=F,quote=F)
save(Untreated_vs_Initial.fc,file="Rdata_files/Untreated_vs_Initial.fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Untreated_vs_Initial.p=data.frame("gene"=row.names(source))
Untreated_vs_Initial.p$lpv=-log(source[,"padj"],10)
Untreated_vs_Initial.p$lpv[source$stat<0]=Untreated_vs_Initial.p$lpv[source$stat<0]*-1
head(Untreated_vs_Initial.p)
write.csv(Untreated_vs_Initial.p,file="Untreated_vs_Initial_lpv.csv",row.names=F,quote=F)
save(Untreated_vs_Initial.p,file="Rdata_files/Untreated_vs_Initial_lpv.RData")


# Treated vs. Initial
# log2 fold changes:
source=Treatment_Treated_vs_Initial[!is.na(Treatment_Treated_vs_Initial$padj),]
Treated_vs_Initial.fc=data.frame("gene"=row.names(source))
Treated_vs_Initial.fc$lfc=source[,"log2FoldChange"]
head(Treated_vs_Initial.fc)
write.csv(Treated_vs_Initial.fc,file="Treated_vs_Initial_fc.csv",row.names=F,quote=F)
save(Treated_vs_Initial.fc,file="Rdata_files/Treated_vs_Initial_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Treated_vs_Initial.p=data.frame("gene"=row.names(source))
Treated_vs_Initial.p$lpv=-log(source[,"padj"],10)
Treated_vs_Initial.p$lpv[source$stat<0]=Treated_vs_Initial.p$lpv[source$stat<0]*-1
head(Treated_vs_Initial.p)
write.csv(Treated_vs_Initial.p,file="Treated_vs_Initial_lpv.csv",row.names=F,quote=F)
save(Treated_vs_Initial.p,file="Rdata_files/Treated_vs_Initial_lpv.RData")


#Treated vs. Untreated
# log2 fold changes:
source=Treatment_Treated_vs_Untreated[!is.na(Treatment_Treated_vs_Untreated$padj),]
Treated_vs_Untreated.fc=data.frame("gene"=row.names(source))
Treated_vs_Untreated.fc$lfc=source[,"log2FoldChange"]
head(Treated_vs_Untreated.fc)
write.csv(Treated_vs_Untreated.fc,file="Treated_vs_Untreated_fc.csv",row.names=F,quote=F)
save(Treated_vs_Untreated.fc,file="Rdata_files/Treated_vs_Untreated_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Treated_vs_Untreated.p=data.frame("gene"=row.names(source))
Treated_vs_Untreated.p$lpv=-log(source[,"padj"],10)
Treated_vs_Untreated.p$lpv[source$stat<0]=Treated_vs_Untreated.p$lpv[source$stat<0]*-1
head(Treated_vs_Untreated.p)
write.csv(Treated_vs_Untreated.p,file="Treated_vs_Untreated_lpv.csv",row.names=F,quote=F)
save(Treated_vs_Untreated.p,file="Rdata_files/Treated_vs_Untreated_lpv.RData")



#### ANNOTATING DGES ####
load("RData_files/realModels_Acer.RData")
load("RData_files/pvals.RData")

library(tidyverse)

#Untreated vs. Initial
Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("UntreatedvsInitial_annotatedDGEs.csv")
#7274 genes 

#Treated vs. Initial
Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("TreatedvsInitial_annotatedDGEs.csv")
#5736 genes


#Treated vs. Untreated
Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("TreatedvsUntreated_annotatedDGEs.csv")
#1588 genes

