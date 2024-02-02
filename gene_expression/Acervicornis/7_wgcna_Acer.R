#### PACKAGES ####

# installing WGCNA:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("flashClust")
# install.packages("WGCNA",dependencies=TRUE)
# repos="http://cran.us.r-project.org"
# run these above commands once, then comment out

# always run these before running any of the following script chunks
library(tidyverse)
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

#### DATA IMPORT and TRAITS ####

# importing data generated from DESeq2 script
lnames=load("RData_files/data4wgcna.RData")
lnames # "vsd.wg"  "design" # log-transformed variance-stabilized gene expression, and table or experimental conditions
datt=t(vsd.wg)
ncol(datt) #15526
nrow(datt) #44

head(design)
str(design)
colnames(vsd.wg)

all.equal(colnames(vsd.wg), rownames(design)) #FALSE

#change rownames to match colnames

rownames(design) <- colnames(vsd.wg)

head(design)

all.equal(colnames(vsd.wg), rownames(design)) #TRUE

rownames(datt)

#change treatment to be binary (0 = FALSE, 1 = TRUE)
Initial = as.numeric(design$Treatment=="Initial")
Untreated = as.numeric(design$Treatment == "Untreated")
Treated = as.numeric(design$Treatment == "Treated")

# assembling table of traits

# coding genotype as binary (0 = FALSE, 1 = TRUE)
# design$Genotype
# SI_C = as.numeric(design$Genotype == "SI_C")
# BC_8b = as.numeric(design$Genotype == "BC_8b")
# MB_B = as.numeric(design$Genotype == "MB_B")

#traits <- data.frame(cbind(Initial, Untreated, Treated, SI_C,BC_8b, MB_B))

traits <- data.frame(cbind(Initial, Untreated, Treated))

rownames(traits) <- rownames(design)

traits

traits %>% 
  rownames_to_column() %>% 
  mutate(rowname = str_extract(rowname, "^[^_]*")) %>% 
  separate(rowname, into = c("Species", "ID"), sep = "\\.") %>% 
  mutate(Species = "Acropora cervicornis") %>% 
  mutate(ID = as.double(ID)) -> traits

traits

days_to_removed <- read_csv("physiotraits_for_WGCNA/days_to_removed.csv")
days_to_removed %>% 
  select(!c("Colony", "Treatment")) %>% 
  filter(Species == "Acropora cervicornis") -> days_to_removed

phagocytosis <- read_csv("physiotraits_for_WGCNA/phagocytosis_alltimepoints.csv")
phagocytosis %>% 
  pivot_wider(names_from="TimePoint", values_from="mean_replicate_percent_perID") %>% 
  dplyr::rename(cells_initial = T0) %>% 
  dplyr::rename(cells_endoftreatment = T2) %>% 
  select(!c("Genotype", "Treatment", "num_days")) %>% 
  filter(Species == "Acer") %>% 
  mutate(Species = "Acropora cervicornis") %>% 
  mutate(ID = gsub("[AP]", "", ID)) %>% 
  mutate(ID = as.double(ID)) -> phagocytosis

Rscore <- read_csv("physiotraits_for_WGCNA/treatment_Rscore.csv")
Rscore %>% 
  filter(Species == "Acropora cervicornis") %>% 
  select(!c("Colony", "Treatment")) %>% 
  mutate(ID = as.double(ID)) -> Rscore

CBASS_fvfm <- read_csv("physiotraits_for_WGCNA/cbass_fvfm_forwgcna.csv")
CBASS_fvfm %>% 
  filter(Species == "Acervicornis") %>% 
  mutate(Species = "Acropora cervicornis") %>% 
  dplyr::rename(ID = Puck) %>% 
  mutate(ID = gsub("[AP]", "", ID)) %>% 
  mutate(ID = as.double(ID)) %>% 
  dplyr::select(!c("Colony", "Treatment")) -> CBASS_fvfm

full_join(traits, days_to_removed, by = c("Species", "ID")) %>% 
  full_join(., phagocytosis) %>% 
  full_join(., Rscore) %>% 
  full_join(., CBASS_fvfm) %>% 
  drop_na(Initial) %>% 
  dplyr::select(!c("Species", "ID")) -> traits_withphysio

traits_withphysio %>% 
  select(!c(1)) -> traits_withphysio

#### OUTLIER DETECTION ####

# identifies outlier genes
gsg = goodSamplesGenes(datt, verbose = 3);
gsg$allOK #if TRUE, no outlier genes
#TRUE!

# calculates mean expression per array, then the number of missing values per array
meanExpressionByArray=apply( datt,1,mean, na.rm=T)
NumberMissingByArray=apply( is.na(data.frame(datt)),1, sum)
NumberMissingByArray
# keep samples with missing values under 500
# in this case, all samples OK

# plots mean expression across all samples
quartz()
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:44), cex.names = 0.7)
# look for any obvious deviations in expression across samples

# sample dendrogram and trait heat map showing outliers
A=adjacency(t(datt),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
#I am going with signed because if you pick unsigned, it mixes negatively and positiviely correlated nodes together and 
#direction of correlation does matter for downstream analysis. 
#(full explanation here: https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/)

# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(traits_withphysio,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(traits_withphysio))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
quartz()
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
# the resulting plot shows a sample dendrogram and the spread of your traits across the cluster
# outlier samples will show as red in the outlierC row

# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datt=datt[!remove.samples,]
traits_withphysio=traits_withphysio[!remove.samples,] #1 sample removed

str(datt) #43 
str(traits_withphysio) #43
  
write.csv(traits, file="traits.csv")

write.csv(traits_withphysio, file="traits_withphysio.csv")

save(datt,traits,file="wgcnaData.RData")


#### SOFT THRESHOLDS ####

library(tidyverse)
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

load("RData_files/wgcnaData.RData")
traits_withphysio <- read_csv("traits_withphysio.csv")

#datt=t(vsd.wg)
str(datt) #43

# Try different betas ("soft threshold") - power factor for calling connections between genes
powers = c(seq(from = 2, to=26, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 8,networkType="signed")

# Plot the results:
# Run from the line below to dev.off()
sizeGrWindow(9, 5)
quartz()
pdf("soft_threshold_signed.pdf",height=4, width=8)

par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#### MAKING MODULES ####

# take a look at the threshold plots produced above, and the output table from the pickSoftThreshold command
# pick the power that corresponds with a SFT.R.sq value above 0.90

# run from the line below to the save command
s.th=19 # re-specify according to previous section
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
rm(adjacency) #for memory space
dissTOM = 1-TOM
rm(TOM)
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average")

plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules.RData")


#### MERGING MODULES ####

mm=load('RData_files/1stPassModules.RData')
mm
lnames=load('RData_files/wgcnaData.RData')
traits
head(datt)

quartz()

MEDissThres = 0.4 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed.RData")


#### MODULE CORRELATIONS ####
# plotting correlations with traits:
load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/wgcnaData.RData");

# Define numbers of genes and samples
nGenes = ncol(datt); #15526
nSamples = nrow(datt); #43
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits_withphysio, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# module-trait correlations

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traits_withphysio),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

table(moduleColors) # gives numbers of genes in each module

# shows only significant correlations
quartz()
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor,2)
textMatrix = cors;
# paste(cors, "\n(",ps, ")", sep = "");
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits_withphysio),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-0.7,0.7),
               main = paste("A.cervicornis Module-Trait correlations"))

# module size barplot
quartz()
labelShift=750 # increase to move module size labels to the right
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,16000),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

# If it was first pass with no module merging, this is where you examine your heatmap and dendrogram of module eigengenes 
#to see where you would like to set cut height (MEDissThres parameter) 
#in the previous section to merge modules that are telling the same story for your trait data 
# A good way to do it is to find a group of similar modules in the heat map and then see at which tree height they connect in the dendrogram

#### GO BACK AND MERGE ####
#done!

#### MODULE MEMBERSHIP SCATTERPLOTS ####

# scatterplots of gene significance (correlation-based) vs kME
load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/wgcnaData.RData");
traits
table(moduleColors)

# run for each of these statements individually
#whichTrait="Initial"
#whichTrait="Treated"
whichTrait="Untreated"

nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(traits_withphysio[,whichTrait]);
names(selTrait) = whichTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");

# selecting specific modules to plot (change depending on which trait you're looking at)
#moduleCols=c("blue","royalblue", "lightgreen", "grey60", "turquoise", "grey") # for Initial
#moduleCols=c("blue","royalblue", "darkgreen", "darkturquoise","grey60", "grey") # for Treated
moduleCols=c("blue","greenyellow","midnightblue", "darkturquoise") # for Untreated

quartz()
# set par to be big enough for all significant module correlations, 
#then run the next whichTrait and moduleCols statements above and repeat from the 'for' loop
par(mfrow=c(1,6)) 

counter=0
# shows correlations for all modules
for(module in modNames[1:length(modNames)]){
counter=counter+1}

quartz()
	par(mfrow=c(3,3))

# shows correlations for significant modules only as specified above
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste(module,"module membership"),
ylab = paste("GS for", whichTrait),
col = "grey50",mgp=c(2.3,1,0))
}


#### EIGENGENE SANITY CHECK ####

# eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules
load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/wgcnaData.RData");

# run for each of these statements individually
#which.module="turquoise"
#which.module="darkturquoise"
#which.module="royalblue"
#which.module="lightgreen"
#which.module="darkgreen"
#which.module="greenyellow"
#which.module="midnightblue"
#which.module="blue"
#which.module="grey"
which.module="grey60"


datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module

#turquoise = 4829
#darkturquoise = 115
#royalblue = 148
#lightgreen = 199
#darkgreen = 118
#greenyellow = 624
#midnightblue = 626
#blue = 4557
#grey = 602
#grey60 = 220


# If individual samples appear to be driving expression of significant modules, they are likely outliers

#### GO/KOG EXPORT ####

# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)
library(WGCNA)
load(file = "RData_files/networkdata_signed.RData") # moduleColors, MEs
load(file = "RData_files/wgcnaData.RData") # vsd table
load(file = "RData_files/data4wgcna.RData") # vsd table

# calculating module memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs)) 
names(allkME)=gsub("kME","",names(allkME))

# run for each of these statements individually

#which.module="turquoise"
#which.module="darkturquoise"
#which.module="royalblue"
#which.module="lightgreen"
#which.module="darkgreen"
#which.module="greenyellow"
#which.module="midnightblue"
#which.module="blue"
#which.module="grey"
which.module="grey60"

# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==which.module)
combo=data.frame("gene"=row.names(t(datt)),"Fish_kME"=allkME[,which.module]*inModuleBinary)
write.csv(combo,file=paste(which.module,".csv",sep=""),row.names=F,quote=F)


#### HEATMAPS ####

# plotting heatmap for named top-kME genes
library(WGCNA)
library(pheatmap)

load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/data4wgcna.RData") 
load(file = "RData_files/wgcnaData.RData");

allkME =signedKME(datt, MEs)
gg=read.delim(file="bioinformatics/Acervicornis_iso2geneName.tab",sep="\t")

#which.module="turquoise"
#which.module="darkturquoise"
#which.module="royalblue"
#which.module="lightgreen"
#which.module="darkgreen"
which.module="greenyellow"
#which.module="midnightblue"
#which.module="blue"
#which.module="grey"
#which.module="grey60"

vsd.wg = t(datt)

top=30 # number of named top-kME genes to plot

datME=MEs
datExpr=datt
modcol=paste("kME",which.module,sep="")
sorted=vsd.wg[order(allkME[,modcol],decreasing=T),]
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
	if (row.names(sorted)[i] %in% gg[,1]) { 
		counts=counts+1
		gn=gg[gg[,1]==row.names(sorted)[i],2]
		gn=paste(gn,row.names(sorted)[i],sep=".")
		if (gn %in% gnames) {
			gn=paste(gn,counts,sep=".")
		}
		gnames=append(gnames,gn) 
		hubs=data.frame(rbind(hubs,sorted[i,]))
		if (counts==top) {break}
	}
} 
row.names(hubs)=gnames

colnames(hubs) 

categories <- c("Initial", "Treated", "Untreated")

# Extract and sort columns for each category
category_columns <- lapply(categories, function(cat) {
  matching_columns <- grep(cat, names(hubs), value = TRUE)
  hubs[, matching_columns, drop = FALSE]
})

# Bind the columns back together in the desired order
reordered_df <- do.call(cbind, category_columns)

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

#pdf(file="heatmap_top30_turquoise.pdf")
pheatmap(as.matrix(reordered_df),scale="row",col=contrasting,border_color=NA,treeheight_col=0,cex=0.9,cluster_rows = F, cluster_cols = F) 
dev.off()


#### HUB GENES ####

library(WGCNA)
library(tidyverse)
load(file = "networkdata_signed.RData")
load(file = "data4wgcna.RData") 
load(file = "wgcnaData.RData");
allkME =as.data.frame(signedKME(datt, MEs))

hubgenes <- chooseTopHubInEachModule(datt, moduleColors, omitColors = "grey", 
                                     power = 2, 
                                     type = "signed")
hubgenes <-data.frame(hubgenes)
hubgenes <- tibble::rownames_to_column(hubgenes, "module")
hubgenes

hubgenes %>%
  rename("gene" = 
           hubgenes) %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> hubgenes
hubgenes

write.csv(hubgenes, file="hubgenes.csv")


