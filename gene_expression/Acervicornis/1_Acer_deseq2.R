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

#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("Rdata_files/vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_fullmodel.pdf", width=15, height=15)
head(vsd)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
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

# plotting PCoA by treatment and Genotype (axes 1 and 2)
pdf(file="plots/PCoA_TreatmentGenotype.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("grey","red","blue")[as.numeric(as.factor(conditions$Treatment))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Genotype)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$Treatment, label=F, col=c("grey","red","blue"))
legend("topright", legend=c("Initial", "Treated", "Untreated"), fill = c("grey","red","blue"), bty="n")
legend("topleft", legend=c("BC-8b", "MB-B", "SI-C"), pch=c(15,17,25), bty="n")
plot(scores[,1], scores[,2],col=c("darkgreen","orange", "black")[as.numeric(as.factor(conditions$Genotype))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Treatment)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Genotype")
ordispider(scores, conditions$Genotype, label=F, col=c("darkgreen","orange", "black"))
legend("topleft", legend=c("BC-8b", "MB-B", "SI-C"), fill = c("darkgreen","orange", "black"), bty="n")
legend("topright", legend=c("Initial", "Treated", "Untreated"), pch=c(15,17,25), bty="n")
dev.off()

# plotting PCoA by treatment
pdf(file="plots/PCoA.pdf")
plot(scores[,1], scores[,2],col=c("grey","blue","red")[as.numeric(as.factor(conditions$Treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$Treatment, label=F, col=c("grey","blue","red"))
legend("topleft", legend=c("Initial", "Untreated", "Treated"), fill = c("grey","blue","red"), bty="n")
dev.off()

# plotting PCoA by treatment and genotype (axes 2 and 3)
pdf(file="PCoA_axes23.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,2], scores[,3],col=c("grey","red","blue")[as.numeric(as.factor(conditions$Treatment))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Genotype)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores[,2:3], conditions$Treatment, label=F, col=c("grey","red","blue"))
legend("topright", legend=c("Initial", "Treated", "Untreated"), fill = c("grey","red","blue"), bty="n")
legend("topleft", legend=c("BC-8b", "MB-B", "SI-C"), pch=c(15,17,25), bty="n")
plot(scores[,2], scores[,3],col=c("darkgreen","orange", "black")[as.numeric(as.factor(conditions$Genotype))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Treatment)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Genotype")
ordispider(scores[,2:3], conditions$Genotype, label=F, col=c("darkgreen","orange", "black"))
legend("topleft", legend=c("BC-8b", "MB-B", "SI-C"), fill = c("darkgreen","orange", "black"), bty="n")
legend("topright", legend=c("Initial", "Treated", "Untreated"), pch=c(15,17,25), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
ad=adonis2(t(vsd)~Genotype + Treatment,data=conditions,method="manhattan",permutations=1e6)
ad
summary(ad)
as.data.frame(ad)

# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:4],labels=row.names(as.data.frame(ad)),col=cols,main="Genotype vs treatment")
dev.off()


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

# treatment
Treatment=results(dds) 
summary(Treatment) 

resultsNames(dds)
#primary ones of interest: "Treatment_Treated_vs_Initial",  "Treatment_Untreated_vs_Initial",  "Treatment_Treated_vs_Untreated"

degs_treatment=row.names(Treatment)[Treatment$padj<0.1 & !(is.na(Treatment$padj))]
#list of gene names that are significant at p < 0.1
length(degs_treatment) #8643 genes

# Treated vs. Untreated
Treatment_Treated_vs_Untreated=results(dds,contrast=c("Treatment","Treated","Untreated"))
summary(Treatment_Treated_vs_Untreated)
degs_Treatment_Treated_vs_Untreated=row.names(Treatment_Treated_vs_Untreated)[Treatment_Treated_vs_Untreated$padj<0.1 & !(is.na(Treatment_Treated_vs_Untreated$padj))]
length(degs_Treatment_Treated_vs_Untreated)

# Treated vs. Initial
Treatment_Treated_vs_Initial=results(dds,contrast=c("Treatment","Treated","Initial"))
summary(Treatment_Treated_vs_Initial)
degs_Treatment_Treated_vs_Initial=row.names(Treatment_Treated_vs_Initial)[Treatment_Treated_vs_Initial$padj<0.1 & !(is.na(Treatment_Treated_vs_Initial$padj))]

# Untreated vs. Initial
Treatment_Untreated_vs_Initial=results(dds,contrast=c("Treatment","Untreated","Initial"))
summary(Treatment_Untreated_vs_Initial)
degs_Treatment_Untreated_vs_Initial=row.names(Treatment_Untreated_vs_Initial)[Treatment_Untreated_vs_Initial$padj<0.1 & !(is.na(Treatment_Untreated_vs_Initial$padj))]

save(Treatment,Treatment_Untreated_vs_Initial,Treatment_Treated_vs_Initial,Treatment_Treated_vs_Untreated, degs_Treatment_Untreated_vs_Initial,degs_Treatment_Treated_vs_Initial, degs_Treatment_Treated_vs_Untreated, file="Rdata_files/pvals.RData")


# density plots: are my DEGs high-abundant or low-abundant?
load("Rdata_files/vsd.RData")
load("Rdata_files/pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density_Treatment.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_Treatment_Untreated_vs_Initial]),col="blue")
lines(density(means[degs_Treatment_Treated_vs_Initial]),col="orange")
lines(density(means[degs_Treatment_Treated_vs_Untreated]),col="lightblue")
legend("topright", title = "Factor", 
       legend=c("Untreated_vs_Initial","Treated_vs_Initial","Treated_vs_Untreated"), 
       fill = c("blue","orange","lightblue"))
dev.off()


#### VENN DIAGRAMS ####

load("RData_files/pvals.RData")
library(DESeq2)

pairwise=list("Untreated vs. Initial"=degs_Treatment_Untreated_vs_Initial, "Treated vs. Initial"=degs_Treatment_Treated_vs_Initial,"Treated vs. Untreated"=degs_Treatment_Treated_vs_Untreated)

# install.packages("VennDiagram")
library(VennDiagram)

# treatment/time contrasts
venn = venn.diagram(
  x = pairwise,  # make sure this contains only three sets
  filename = NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582"),  # three colors, one for each set
  alpha = 0.5,
  label.col = c("red3", "white", "cornflowerblue", "black", "white", "white", "black"),  # colors for each of the 7 regions
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "red3"),  # three category colors
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0, 0.5), c(0.75, 0.5), c(0.5, 0.5))  # adjust positions for three categories
)

           
pdf(file="Venn_Acer.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


library(ggvenn)

ggvenn(pairwise) + 
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

#Untreated vs. Initial
as.data.frame(Treatment_Untreated_vs_Initial) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("UntreatedvsInitial_annotatedDGEs.csv")
#8,643 genes 

#Treated vs. Initial
as.data.frame(Treatment_Treated_vs_Initial) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("TreatedvsInitial_annotatedDGEs.csv")
#7,004 genes


#Treated vs. Untreated
as.data.frame(Treatment_Treated_vs_Untreated) %>%
  rownames_to_column(var="gene") %>% 
  filter(padj < 0.1) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("TreatedvsUntreated_annotatedDGEs.csv")
#2,295 genes


#### CHERRY PICKING ####

load("RData_files/Untreated_vs_Initial_lpv.RData")
Untreated_vs_Initial.p %>% 
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% 
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) %>% 
  write.csv("UntreatedvsInitial_cherrypicking.csv")

load("RData_files/Treated_vs_Initial_lpv.RData")
Treated_vs_Initial.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% 
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) %>% 
  write.csv("TreatedvsInitial_cherrypicking.csv")

load("RData_files/Treated_vs_Untreated_lpv.RData")
Treated_vs_Untreated.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% 
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) %>% 
  write.csv("TreatedvsUntreated_cherrypicking.csv")



