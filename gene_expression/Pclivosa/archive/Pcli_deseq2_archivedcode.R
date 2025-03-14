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