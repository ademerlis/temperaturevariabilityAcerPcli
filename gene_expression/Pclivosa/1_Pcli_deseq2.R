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
#install.packages("ggvenn")
#install.packages("ggforce")


#### DATA IMPORT ####
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)
library(tidyverse)

#read in counts
counts = readxl::read_xlsx("../MS_Bioinformatics/MS_bioinformatics_alignmenttests/pcli/magana/allcounts_pcli.xlsx") #alignment to Avila-Magana et al. 2021 transcriptome 

column_to_rownames(counts, var ="...1") -> counts

colnames(counts) # sample IDs (i.e. Pcli-003). You need to keep track of column names and row names of data frames throughout this analysis because they MUST MATCH 
rownames(counts) #gene names (but not actual gene names, arbitrary gene names. To get "actual" gene names, need annotation files)

# how many genes we have total?
nrow(counts) #59947
ncol(counts) #48 samples

# how does the data look? 
head(counts)

# filtering out low-count genes
keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData) #34300
ncol(countData) #48
#write.csv(countData, file = "results_csv/Pcli_countdata.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna) #13867
ncol(counts4wgcna) #48
#write.csv(counts4wgcna, file="results_csv/Pcli_counts4wgcna.csv")

# importing a design .csv file
design = read.csv("../../treatment_metadata.csv", head=TRUE)

design %>% 
  select(Species:Treatment) %>% 
  filter(Species == "Pcli") -> design

design$Genotype <- as.factor(design$Genotype)
design$Treatment <- as.factor(design$Treatment)
column_to_rownames(design, var="Sample_ID") -> design
str(design)


#### FULL MODEL DESIGN (Genotype + Treatment) and OUTLIERS ####

#when making dds formula, it is CRITICAL that you put the right order of variables in the design. The design indicates how to model the samples, 
#(here: design = ~batch + condition), 
#that we want to measure the effect of the condition, controlling for batch differences. The two factor variables batch and condition should be columns of coldata.

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)

# reorders fate factor according to "control" vs "treatment" levels
dds$Treatment <- factor(dds$Treatment, levels = c("Initial", "Untreated", "Treated"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

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
outs=c(5,23,39,46) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)
dds$Treatment <- factor(dds$Treatment, levels = c("Initial", "Untreated", "Treated"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="Rdata_files/initial_fullddsdesigncountsVsdcountsWGCNA.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
colnames(vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2],design[,5],sep=".")
snames
# rename the column names
colnames(vsd)=snames

save(vsd,design,file="Rdata_files/vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ Genotype + Treatment)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=FALSE)  #blind=TRUE is the default, and it is a fully unsupervised transformation. However, the creator of DESeq2,
#Michael Love, recommends using blind=FALSE for downstream analyses because when transforming data, the full use of the design information should be made. If many genes have
#large differences in counts due to experimental design, then blind=FALSE will account for that.

head(vsd.wg)
colnames(vsd.wg)=snames
save(vsd.wg,design,file="Rdata_files/data4wgcna.RData")


#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("Rdata_files/vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="plots/heatmap_fullmodel.pdf", width=15, height=15)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
conditions$Treatment <- factor(conditions$Treatment, levels = c("Initial", "Untreated", "Treated"))

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
pdf(file="plots/PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's
#there is 1 good PC based on this figure

# plotting PCoA by treatment
pdf(file="plots/PCoA.pdf")
plot(scores[,1], scores[,2],col=c("black","blue","red")[as.numeric(as.factor(conditions$Treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$Treatment, label=F, col=c("black","blue","red"))
legend("topleft", legend=c("Initial", "Untreated", "Treated"), fill = c("black","blue","red"), bty="n")
dev.off()

# plotting PCoA by treatment and Genotype (axes 1 and 2)
pdf(file="plots/PCoA_TreatmentGenotype.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("grey","red","blue")[as.numeric(as.factor(conditions$Treatment))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Genotype)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$Treatment, label=F, col=c("grey","red","blue"))
legend("topright", legend=c("Initial", "Treated", "Untreated"), fill = c("grey","red","blue"), bty="n")
legend("topleft", legend=c("A", "B", "C"), pch=c(15,17,25), bty="n")
plot(scores[,1], scores[,2],col=c("darkgreen","orange", "black")[as.numeric(as.factor(conditions$Genotype))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Treatment)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Genotype")
ordispider(scores, conditions$Genotype, label=F, col=c("darkgreen","orange", "black"))
legend("topleft", legend=c("A", "B", "C"), fill = c("darkgreen","orange", "black"), bty="n")
legend("topright", legend=c("Initial", "Treated", "Untreated"), pch=c(15,17,25), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="plots/PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matrices: 
ad=adonis2(t(vsd)~Genotype + Treatment,data=conditions,method="manhattan",permutations=1e6)
ad
summary(ad)
as.data.frame(ad)

# creating pie chart to represent ANOVA results
cols=c("orange", "blue","lightblue","grey80")
pdf(file="plots/ANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:4],labels=row.names(as.data.frame(ad)),col=cols,main="Genotype vs Treatment")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design
load("initial_fullddsdesigncountsVsdcountsWGCNA.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements

#dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Treatment)
rownames(design)
colnames(countData)
dds=DESeq(dds, parallel=TRUE)

# saving all models
save(dds,file="Rdata_files/realModels_Pcli.RData")


#### DEGs and CONTRASTS ####

load("Rdata_files/realModels_Pcli.RData")
library(DESeq2)
library(tidyverse)

read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
           sep = "\t",
           quote="", fill=FALSE) %>%
  rename(gene = V1,
         annot = V2) -> iso2geneName

# treatment
Treatment=results(dds) 
summary(Treatment, alpha = 0.05) 
degs_Treatment=row.names(Treatment)[Treatment$padj<0.05 & !(is.na(Treatment$padj))]
length(degs_Treatment) #141

# export for GO analysis
Treatment %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene, padj) %>% 
  write_csv("degs_Treatment.csv")

resultsNames(dds)
#"Intercept"                      "Genotype_B_vs_A"                "Genotype_C_vs_A"               
#"Treatment_Untreated_vs_Initial" "Treatment_Treated_vs_Initial" 

# treatment and time contrasts
Treatment_Untreated_vs_Initial=results(dds,contrast=c("Treatment","Untreated","Initial"))
summary(Treatment_Untreated_vs_Initial, alpha = 0.05)
degs_Treatment_Untreated_vs_Initial=row.names(Treatment_Untreated_vs_Initial)[Treatment_Untreated_vs_Initial$padj<0.05 & !(is.na(Treatment_Untreated_vs_Initial$padj))]
length(degs_Treatment_Untreated_vs_Initial) #132 genes

Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("UntreatedvInitial_annotDGEs.csv")

# export for GO analysis
Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene, padj) %>% 
  write_csv("degs_Treatment_Untreated_vs_Initial.csv")

Treatment_Treated_vs_Initial=results(dds,contrast=c("Treatment","Treated","Initial"))
summary(Treatment_Treated_vs_Initial, alpha = 0.05)
degs_Treatment_Treated_vs_Initial=row.names(Treatment_Treated_vs_Initial)[Treatment_Treated_vs_Initial$padj<0.05 & !(is.na(Treatment_Treated_vs_Initial$padj))]
length(degs_Treatment_Treated_vs_Initial) #141 genes

Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("TreatedvInitial_annotDGEs.csv")

# export for GO analysis
Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene, padj) %>% 
  write_csv("degs_Treatment_Treated_vs_Initial.csv")

Treatment_Treated_vs_Untreated=results(dds,contrast=c("Treatment","Treated","Untreated"))
summary(Treatment_Treated_vs_Untreated, alpha = 0.05)
degs_Treatment_Treated_vs_Untreated=row.names(Treatment_Treated_vs_Untreated)[Treatment_Treated_vs_Untreated$padj<0.05 & !(is.na(Treatment_Treated_vs_Untreated$padj))]
length(degs_Treatment_Treated_vs_Untreated) #21 genes

Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("TreatedvUntreated_annotDGEs.csv")

# export for GO analysis
Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::select(gene, padj) %>% 
  write_csv("degs_Treatment_Treated_vs_Untreated.csv")


save(Treatment,Treatment_Untreated_vs_Initial,Treatment_Treated_vs_Initial,Treatment_Treated_vs_Untreated, degs_Treatment_Untreated_vs_Initial,degs_Treatment_Treated_vs_Initial, degs_Treatment_Treated_vs_Untreated, file="Rdata_files/pvals.RData")

# density plots: are my DEGs high-abundant or low-abundant?
load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="plots/DEG_density_Treatment.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_Treatment_Untreated_vs_Initial]),col="blue")
lines(density(means[degs_Treatment_Treated_vs_Initial]),col="lightblue")
lines(density(means[degs_Treatment_Treated_vs_Untreated]),col="yellow")
legend("topright", title = "Treatment", 
       legend=c("Untreated_vs_Initial","Treated_vs_Initial","Treated_vs_Untreated"), 
       fill = c("blue","lightblue","yellow"))
dev.off()


#### VENN DIAGRAMS ####

load("Rdata_files/pvals.RData")
library(DESeq2)

pairwise=list("Untreated_vs_Initial"=degs_Treatment_Untreated_vs_Initial, "Treated_vs_Initial"=degs_Treatment_Treated_vs_Initial,"Treated_vs_Untreated"=degs_Treatment_Treated_vs_Untreated)

# Function to find common elements in a list of vectors
find_common_elements <- function(lst) {
  common_elements <- lst[[1]]
  for (vec in lst[-1]) {
    common_elements <- intersect(common_elements, vec)
  }
  return(common_elements)
}

# Find common elements
common_elements <- find_common_elements(pairwise)

# Print the common elements
print(common_elements) #there are no shared genes

library(ggvenn)

ggvenn(pairwise) + 
  scale_fill_manual(values = c("#ca0020", "#0571b0", "#f4a582"))


library(VennDiagram)

# treatment contrasts
venn = venn.diagram(
  x = pairwise,  # Make sure this is now for three groups
  filename = NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582"),  # Adjusted for three groups
  alpha = 0.5,
  label.col = c("red3", "cornflowerblue", "black", "white", "black", "darkred", "black"),  # Adjusted for three groups
  fontface = "bold",
  cat.col = c("darkred", "darkblue", "red3"),  # Adjusted for three groups
  cex = 2.0,  # Adjusts the size of the text in the diagram, including region labels
  cat.default.pos = "text",  # Set category labels to be outside the circles
  cat.cex = 2.5,  # Adjusts the size of the category labels
  scale = 0.5,  # Adjust this value to scale the diagram
  margins = c(0.1, 0.1, 0.1, 0.1) 
  )


pdf(file="Venn_Pcli.pdf", height=10, width=10)
grid.draw(venn)
dev.off()


#### KOG EXPORT ####

load("Rdata_files/realModels_Pcli.RData")
load("Rdata_files/pvals.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# Untreated vs Initial
# log2 fold changes:
source=Treatment_Untreated_vs_Initial[!is.na(Treatment_Untreated_vs_Initial$padj),]
Untreated_vs_Initial.fc=data.frame("gene"=row.names(source))
Untreated_vs_Initial.fc$lfc=source[,"log2FoldChange"]
head(Untreated_vs_Initial.fc)
write.csv(Untreated_vs_Initial.fc,file="results_csv/Untreated_vs_Initial_fc.csv",row.names=F,quote=F)
save(Untreated_vs_Initial.fc,file="Rdata_files/Untreated_vs_Initial.fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Untreated_vs_Initial.p=data.frame("gene"=row.names(source))
Untreated_vs_Initial.p$lpv=-log(source[,"padj"],10)
Untreated_vs_Initial.p$lpv[source$stat<0]=Untreated_vs_Initial.p$lpv[source$stat<0]*-1
head(Untreated_vs_Initial.p)
write.csv(Untreated_vs_Initial.p,file="results_csv/Untreated_vs_Initial_lpv.csv",row.names=F,quote=F)
save(Untreated_vs_Initial.p,file="Rdata_files/Untreated_vs_Initial_lpv.RData")


# Treated vs. Initial
# log2 fold changes:
source=Treatment_Treated_vs_Initial[!is.na(Treatment_Treated_vs_Initial$padj),]
Treated_vs_Initial.fc=data.frame("gene"=row.names(source))
Treated_vs_Initial.fc$lfc=source[,"log2FoldChange"]
head(Treated_vs_Initial.fc)
write.csv(Treated_vs_Initial.fc,file="results_csv/Treated_vs_Initial_fc.csv",row.names=F,quote=F)
save(Treated_vs_Initial.fc,file="Rdata_files/Treated_vs_Initial_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Treated_vs_Initial.p=data.frame("gene"=row.names(source))
Treated_vs_Initial.p$lpv=-log(source[,"padj"],10)
Treated_vs_Initial.p$lpv[source$stat<0]=Treated_vs_Initial.p$lpv[source$stat<0]*-1
head(Treated_vs_Initial.p)
write.csv(Treated_vs_Initial.p,file="results_csv/Treated_vs_Initial_lpv.csv",row.names=F,quote=F)
save(Treated_vs_Initial.p,file="Rdata_files/Treated_vs_Initial_lpv.RData")


#Treated vs. Untreated
# log2 fold changes:
source=Treatment_Treated_vs_Untreated[!is.na(Treatment_Treated_vs_Untreated$padj),]
Treated_vs_Untreated.fc=data.frame("gene"=row.names(source))
Treated_vs_Untreated.fc$lfc=source[,"log2FoldChange"]
head(Treated_vs_Untreated.fc)
write.csv(Treated_vs_Untreated.fc,file="results_csv/Treated_vs_Untreated_fc.csv",row.names=F,quote=F)
save(Treated_vs_Untreated.fc,file="Rdata_files/Treated_vs_Untreated_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
Treated_vs_Untreated.p=data.frame("gene"=row.names(source))
Treated_vs_Untreated.p$lpv=-log(source[,"padj"],10)
Treated_vs_Untreated.p$lpv[source$stat<0]=Treated_vs_Untreated.p$lpv[source$stat<0]*-1
head(Treated_vs_Untreated.p)
write.csv(Treated_vs_Untreated.p,file="results_csv/Treated_vs_Untreated_lpv.csv",row.names=F,quote=F)
save(Treated_vs_Untreated.p,file="Rdata_files/Treated_vs_Untreated_lpv.RData")


#### ANNOTATING DGES ####

#using the log-10 transformed p-adj value because that's what's used in the 6_commongenes.R script and I want to make sure all the numbers match

#Untreated vs. Initial
Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("Untreated_vs_Initial_annotatedDGEs_05.csv")
#132 genes 

# Treated vs. Initial
Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("Treated_vs_Initial_annotatedDGEs_05.csv")
#141 genes


# Treated vs. Untreated
Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("Treated_vs_Untreated_annotatedDGEs_05.csv")
#21 genes

