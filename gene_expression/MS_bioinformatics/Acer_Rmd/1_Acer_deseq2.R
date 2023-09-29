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
counts = readxl::read_xlsx("../Bioinformatics/acer/libro/allcounts_acer.xlsx") #Libro et al. 2013 genome had the best alignment rate, so that is what we will use

column_to_rownames(counts, var ="...1") -> counts

# how many genes we have total?
nrow(counts) #57358
ncol(counts) #48 samples

# how does the data look? 
head(counts)

colnames(counts) # sample IDs (i.e. Acer-005). You need to keep track of column names and row names of data frames throughout this analysis because they MUST MATCH 
rownames(counts) #gene names (but not actual gene names, arbitrary gene names like Acropora000001. To get "actual" gene names, need annotation files)

# filtering out low-count genes
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData) #47882
ncol(countData) #48
#write.csv(countData, file = "Acer_countdata.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna) #27829
ncol(counts4wgcna) #48
#write.csv(counts4wgcna, file="Acer_counts4wgcna.csv")

colnames(counts4wgcna) #still same sample ID names "Acer-005" etc

# importing a design .csv file
design = read.csv("../../RNA_extraction_sequencing_data.csv", head=TRUE)

design %>% 
  select(Species:Treatment) %>% 
  filter(Species == "Acer") -> design

design %>% 
  mutate(time_point = case_when(Experiment.phase == "Pre-treatment" ~ "Day_0",
                                Experiment.phase == "last day of treatment" ~ "Day_29")) -> design
column_to_rownames(design, var="Sample_ID") -> design
design$Genotype <- as.factor(design$Genotype)
design$Genotype <- factor(gsub("-", "_", design$Genotype)) #DESeq2 does not like hyphens in factor names
design$Treatment <- as.factor(design$Treatment)
design$group <- factor(paste0(design$Treatment, "_", design$time_point))


colnames(design) #"Species"   "Genotype"   "Sample.ID"  "Date.Sampled" "Experiment.phase" "Treatment"   "time_point"  "group" 
rownames(design) #"Acer-005" etc

#### FULL MODEL DESIGN (Genotype + group) and OUTLIERS ####

#when making dds formula, it is CRITICAL that you put the right order of variables in the design. The design indicates how to model the samples, 
#(here: design = ~batch + condition), 
#that we want to measure the effect of the condition, controlling for batch differences. The two factor variables batch and condition should be columns of coldata.

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + group)

# reorders factor according to "control" vs "treatment" levels
dds$group <- factor(dds$group, levels = c("control_Day_0","control_Day_29","variable_Day_0","variable_Day_29"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

colnames(Vsd)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("group"),force=T) #Genotype is not included as an intgroup because it is not the main factors of interest.
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested.
#I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. 
#Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(20,22,28) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + group)
dds$group <- factor(dds$group, levels = c("control_Day_0","control_Day_29","variable_Day_0","variable_Day_29"))

#These values should all be the same:
colnames(dds)
rownames(design)
colnames(countData)
colnames(counts4wgcna)

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="initial_fullddsdesigncountsVsdcountsWGCNA.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
colnames(vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2], design[,8],sep=".")
snames #i.e. Acer-005.SI_C.control_Day_0

# renames the column names
colnames(vsd)=snames

save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ Genotype + group)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=FALSE) #blind=TRUE is the default, and it is a fully unsupervised transformation. However, the creator of DESeq2,
#Michael Love, recommends using blind=FALSE for downstream analyses because when transforming data, the full use of the design information should be made. If many genes have
#large differences in counts due to experimental design, then blind=FALSE will account for that.

head(vsd.wg)
colnames(vsd.wg)=snames
colnames(vsd.wg)
colnames(vsd)
save(vsd.wg,design,file="data4wgcna.RData")

#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("~/OneDrive - University of Miami/GitHub/Ch2_temperaturevariability2023/gene_expression/MS_bioinformatics/Acer_Rmd/RData_files/vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_fullmodel.pdf", width=15, height=15)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
conditions$group <- factor(conditions$group, levels = c("control_Day_0","control_Day_29","variable_Day_0","variable_Day_29"))

colnames(vsd)
rownames(conditions)
#while the names aren't equal, they are in the same order, so when doing the PCoA and using conditions it still should work

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
# copy this table for % variation explained by each axis (Relative_eig column)
dds.pcoa$values

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's
#there are 3 "good PCs" based on this figure

# plotting PCoA by treatment and time (axes 1 and 2)
pdf(file="PCoA.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("red","blue")[as.numeric(as.factor(conditions$Treatment))],pch=c(1,19)[as.numeric(as.factor(conditions$time_point))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$Treatment, label=F, col=c("red","blue"))
legend("topright", legend=c("Control", "Variable"), fill = c("red","blue"), bty="n")
legend("topleft", legend=c("Day_0", "Day_29"), pch=c(1,19), bty="n")
plot(scores[,1], scores[,2],col=c("grey","black")[as.numeric(as.factor(conditions$time_point))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Treatment)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
ordispider(scores, conditions$time, label=F, col=c("grey","black"))
legend("topleft", legend=c("Day_0", "Day_29"), fill = c("grey","black"), bty="n")
legend("topright", legend=c("Control", "Variable"), pch=c(15,17,25), bty="n")
dev.off()

# plotting PCoA by treatment and time (axes 2 and 3)
pdf(file="PCoA_axes23.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,2], scores[,3],col=c("red","blue")[as.numeric(as.factor(conditions$Treatment))],pch=c(1,19)[as.numeric(as.factor(conditions$time_point))], xlab="Coordinate 2", ylab="Coordinate 3", main="Treatment")
ordispider(scores[,2:3], conditions$Treatment, label=F, col=c("red","blue"))
legend("topright", legend=c("Control", "Variable"), fill = c("red","blue"), bty="n")
legend("topleft", legend=c("Day_0", "Day_29"), pch=c(1,19), bty="n")
plot(scores[,2], scores[,3],col=c("grey","black")[as.numeric(as.factor(conditions$time_point))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Treatment)))], xlab="Coordinate 2", ylab="Coordinate 3", main="Time")
ordispider(scores[,2:3], conditions$time, label=F, col=c("grey","black"))
legend("topleft", legend=c("Day_0", "Day_29"), fill = c("grey","black"), bty="n")
legend("topright", legend=c("Control", "Variable"), pch=c(15,17,25), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
#do the results of the PERMANOVA change depending on whether the column and row names match?
rownames(conditions) = snames 
rownames(conditions)
conditions
ad=adonis2(t(vsd)~Genotype + time_point*Treatment,data=conditions,method="manhattan",permutations=1e6)
ad
summary(ad)
as.data.frame(ad)

# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:4],labels=row.names(as.data.frame(ad)),col=cols,main="time vs treatment")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design
load("initial_fullddsdesigncountsVsdcountsWGCNA.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
#dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + group)
rownames(design)
colnames(countData)
dds=DESeq(dds, parallel=TRUE)

# saving all models
save(dds,file="realModels_Acer.RData")


#### DEGs and CONTRASTS ####

load("realModels_Acer.RData")
library(DESeq2)

# treatment
treatment_time=results(dds) 
summary(treatment_time) 
degs_treatment_time=row.names(treatment_time)[treatment_time$padj<0.1 & !(is.na(treatment_time$padj))]
#list of gene names that are significant at p < 0.1
length(degs_treatment_time) #12,836 genes

resultsNames(dds)
#"Intercept"                              "Genotype_MB_B_vs_BC_8b"                 "Genotype_SI_C_vs_BC_8b"                
#"group_control_Day_29_vs_control_Day_0"  "group_variable_Day_0_vs_control_Day_0"  "group_variable_Day_29_vs_control_Day_0"
#there are more contrasts though, so you have to specify 

# treatment and time contrasts
variable29_control29=results(dds,contrast=c("group","variable_Day_29","control_Day_29"))
summary(variable29_control29)
degs_variable29_control29=row.names(variable29_control29)[variable29_control29$padj<0.1 & !(is.na(variable29_control29$padj))]
length(degs_variable29_control29) #4,060 genes

variable0_control0=results(dds,contrast=c("group","variable_Day_0","control_Day_0"))
summary(variable0_control0)
degs_variable0_control0=row.names(variable0_control0)[variable0_control0$padj<0.1 & !(is.na(variable0_control0$padj))]
length(degs_variable0_control0) #149 genes

variable0_variable29=results(dds,contrast=c("group","variable_Day_0","variable_Day_29"))
summary(variable0_variable29)
degs_variable0_variable29=row.names(variable0_variable29)[variable0_variable29$padj<0.1 & !(is.na(variable0_variable29$padj))]
length(degs_variable0_variable29) #10,540 genes

control0_control29=results(dds,contrast=c("group","control_Day_0","control_Day_29"))
summary(control0_control29)
degs_control0_control29=row.names(control0_control29)[control0_control29$padj<0.1 & !(is.na(control0_control29$padj))]
length(degs_control0_control29) #16,022 genes

save(treatment_time,variable29_control29,variable0_control0,variable0_variable29,control0_control29, degs_variable29_control29,degs_variable0_control0, degs_variable0_variable29, degs_control0_control29, file="pvals.RData")


# density plots: are my DEGs high-abundant or low-abundant?
load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density_treatment.time.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_variable29_control29]),col="blue")
lines(density(means[degs_variable0_control0]),col="orange")
lines(density(means[degs_variable0_variable29]),col="lightblue")
lines(density(means[degs_control0_control29]),col="yellow")
legend("topright", title = "Factor", 
       legend=c("variable29_control29","variable0_control0","variable0_variable29","control0_control29"), 
       fill = c("blue","orange","lightblue","yellow"))
dev.off()


#### VENN DIAGRAMS ####

load("RData_files/pvals.RData")
library(DESeq2)

pairwise=list("V29/C29"=degs_variable29_control29, "V0/C0"=degs_variable0_control0,"V0/V29"=degs_variable0_variable29,"C0/C29"=degs_control0_control29)

# install.packages("VennDiagram")
library(VennDiagram)

# treatment/time contrasts
venn=venn.diagram(
  x = pairwise,
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","black","white","white","black","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
           
pdf(file="Venn_Acer.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


#### GO/KOG EXPORT ####

load("RData_files/realModels_Acer.RData")
load("RData_files/pvals.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# variable 29 vs control 29
# log2 fold changes:
source=variable29_control29[!is.na(variable29_control29$padj),]
variable29_control29.fc=data.frame("gene"=row.names(source))
variable29_control29.fc$lfc=source[,"log2FoldChange"]
head(variable29_control29.fc)
write.csv(variable29_control29.fc,file="variable29_control29_fc.csv",row.names=F,quote=F)
save(variable29_control29.fc,file="variable29_control29_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
variable29_control29.p=data.frame("gene"=row.names(source))
variable29_control29.p$lpv=-log(source[,"padj"],10)
variable29_control29.p$lpv[source$stat<0]=variable29_control29.p$lpv[source$stat<0]*-1
head(variable29_control29.p)
write.csv(variable29_control29.p,file="variable29_control29_lpv.csv",row.names=F,quote=F)
save(variable29_control29.p,file="variable29_control29_lpv.RData")

#length of both variable29_control29 data frames = 30,233 genes. This is not a list of significant genes, 
#but a list of all genes with p-adj and log2foldchange values that are not NA

# variable0 vs control0
# log2 fold changes:
source=variable0_control0[!is.na(variable0_control0$padj),]
variable0_control0.fc=data.frame("gene"=row.names(source))
variable0_control0.fc$lfc=source[,"log2FoldChange"]
head(variable0_control0.fc)
write.csv(variable0_control0.fc,file="variable0_control0_fc.csv",row.names=F,quote=F)
save(variable0_control0.fc,file="variable0_control0_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
variable0_control0.p=data.frame("gene"=row.names(source))
variable0_control0.p$lpv=-log(source[,"padj"],10)
variable0_control0.p$lpv[source$stat<0]=variable0_control0.p$lpv[source$stat<0]*-1
head(variable0_control0.p)
write.csv(variable0_control0.p,file="variable0_control0_lpv.csv",row.names=F,quote=F)
save(variable0_control0.p,file="variable0_control0_lpv.RData")

#length of both variable0_control0 data frames = 20,952 genes. This is not a list of significant genes, 
#but a list of all genes with p-adj and log2foldchange values that are not NA

#variable0 vs variable29
# log2 fold changes:
source=variable0_variable29[!is.na(variable0_variable29$padj),]
variable0_variable29.fc=data.frame("gene"=row.names(source))
variable0_variable29.fc$lfc=source[,"log2FoldChange"]
head(variable0_variable29.fc)
write.csv(variable0_variable29.fc,file="variable0_variable29_fc.csv",row.names=F,quote=F)
save(variable0_variable29.fc,file="variable0_variable29_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
variable0_variable29.p=data.frame("gene"=row.names(source))
variable0_variable29.p$lpv=-log(source[,"padj"],10)
variable0_variable29.p$lpv[source$stat<0]=variable0_variable29.p$lpv[source$stat<0]*-1
head(variable0_variable29.p)
nrows(variable0_variable29.p)
write.csv(variable0_variable29.p,file="variable0_variable29_lpv.csv",row.names=F,quote=F)
save(variable0_variable29.p,file="variable0_variable29_lpv.RData")

#length of both variable0_variable29 data frames = 36,731 genes. This is not a list of significant genes, 
#but a list of all genes with p-adj and log2foldchange values that are not NA

# control0_control29
# log2 fold changes:
source=control0_control29[!is.na(control0_control29$padj),]
control0_control29.fc=data.frame("gene"=row.names(source))
control0_control29.fc$lfc=source[,"log2FoldChange"]
head(control0_control29.fc)
write.csv(control0_control29.fc,file="control0_control29_fc.csv",row.names=F,quote=F)
save(control0_control29.fc,file="control0_control29_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
control0_control29.p=data.frame("gene"=row.names(source))
control0_control29.p$lpv=-log(source[,"padj"],10)
control0_control29.p$lpv[source$stat<0]=control0_control29.p$lpv[source$stat<0]*-1
head(control0_control29.p)
write.csv(control0_control29.p,file="control0_control29_lpv.csv",row.names=F,quote=F)
save(control0_control29.p,file="control0_control29_lpv.RData")

#length of both control0_control29 data frames = 40,443 genes. This is not a list of significant genes, 
#but a list of all genes with p-adj and log2foldchange values that are not NA


#### ANNOTATING DGES ####
load("RData_files/realModels_Acer.RData")
load("RData_files/pvals.RData")

#Control_Day0 vs Control_Day29
as.data.frame(control0_control29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("control0_control29_annotatedDGEs.csv")
#16,022 genes 

as.data.frame(control0_control29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1 & abs(log2FoldChange) >2) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("control0_control29_annot_padjL2FC.csv")
#864 genes

#Variable Day 0 vs Variable Day 29
as.data.frame(variable0_variable29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable0_variable29_annotatedDGEs.csv")
#10,540 genes

as.data.frame(variable0_variable29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1 & abs(log2FoldChange) >= 2) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable0_variable29_annot_padjL2FC.csv")
#273 genes

#Variable_Day29 vs Control_Day29
as.data.frame(variable29_control29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable29_control29_annotatedDGEs.csv")
#4,060 genes

as.data.frame(variable29_control29) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1 & abs(log2FoldChange) >= 2) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable29_control29_annot_padjL2FC.csv")
#126 genes

#Variable Day 0 vs Control Day 0
as.data.frame(variable0_control0) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable0_control0_annotatedDGEs.csv")
#149 genes

as.data.frame(variable0_control0) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1 & abs(log2FoldChange) >= 2) %>% 
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("variable0_control0_annot_padjL2FC.csv")
# 1 gene


#### CHERRY PICKING ####

load("RData_files/variable0_control0_lpv.RData")
variable0_variable29.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% 
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) %>% 
  write.csv("variable0variable29_cherrypicking.csv")

load("RData_files/variable29_control29_lpv.RData")
control0_control29.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% 
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) %>% 
  write.csv("control0control29_cherrypicking.csv")



