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
ad=adonis2(t(vsd2)~Genotype + Treatment,data=design,method="manhattan",permutations=1e6)
ad
summary(ad)
as.data.frame(ad)

# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:4],labels=row.names(as.data.frame(ad)),col=cols,main="Genotype vs treatment")
dev.off()

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




