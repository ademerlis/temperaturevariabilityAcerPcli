# Title: A. cervicornis GFF adjustments
# Project: Sedimentation RNA-Seq / Mcap 2020
# Author: J. Ashey / A. Huffmyer -> Allyson DeMerlis
# Date: 06/28/2023

# Need to do some acerv gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.table(file="~/Downloads/genomic.gff", header=FALSE, sep="\t") 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 
Acerv.gff$transcript_id <- gsub("Parent=", "", Acerv.gff$transcript_id) #remove Parent=
head(Acerv.gff)

# Create Parent ID 
Acerv.gff$parent_id <- sub(".*Parent=", "", Acerv.gff$gene)
Acerv.gff$parent_id <- sub(";.*", "", Acerv.gff$parent_id)
Acerv.gff$parent_id <- gsub("ID=", "", Acerv.gff$parent_id) #remove ID= 
head(Acerv.gff)

Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id, ";gene_id=", Acerv.gff$parent_id),  paste0(gene)))
head(Acerv.gff)

Acerv.gff<-Acerv.gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

#save file
write.table(Acerv.gff, file="~/Acer_2023.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

