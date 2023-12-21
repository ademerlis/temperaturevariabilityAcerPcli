#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### DESEQ IMPORT ####
load("RData_files/pvals.RData")

load("RData_files/Treated_vs_Initial_lpv.RData")
Treated_vs_Initial.p <- Treated_vs_Initial.p %>% dplyr::rename("lpv_Treated_vs_Initial" = lpv)
str(Treated_vs_Initial.p) #16,333 genes

load("RData_files/Untreated_vs_Initial_lpv.RData")
Untreated_vs_Initial.p <- Untreated_vs_Initial.p %>% dplyr::rename("lpv_Untreated_vs_Initial" = lpv)
str(Untreated_vs_Initial.p) #17,662 genes

load("RData_files/Treated_vs_Untreated_lpv.RData")
Treated_vs_Untreated.p <- Treated_vs_Untreated.p %>% dplyr::rename("lpv_Treated_vs_Untreated" = lpv) 
str(Treated_vs_Untreated.p) #17,662 genes


#### DEG MATCHING ####

# These sections of code do several things: 1) join common DEGs across experiments with -log10(padj), 
#2) filter by 0.1 padj cutoff (log10(0.1)=1), 3)  adds gene annotations, and 4) then pulls on corresponding KOG classes

# NOTE: normally a -log10 transformation of a p-value between 0 and 1 would not result in a negative number. However, when the above files 
#were generated, (i.e. control0_control29.p), the stat of the directionality of the expression pattern (i.e. negative or positive) was encoded
#in the lpv value manually. The exact lines of code for the transformations were:
#variable29_control29.p$lpv=-log(source[,"padj"],10)
#variable29_control29.p$lpv[source$stat<0]=variable29_control29.p$lpv[source$stat<0]*-1
#This is why the absolute value of lpv is taken for obtaining genes of significance >= 1 (which corresponds to alpha = 0.1)


# first see if there are any shared genes with Treated_vs_Initial and Untreated_vs_Initial
Treated_vs_Initial.p %>%
  inner_join(Untreated_vs_Initial.p, by = "gene") %>%
  filter(abs(lpv_Treated_vs_Initial) >= 1 & abs(lpv_Untreated_vs_Initial) >= 1) %>%
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> TreatedvInitial_vs_UntreatedvInitial 

str(TreatedvInitial_vs_UntreatedvInitial)
#there are 58 shared genes (cut-off is p-adj of 0.1)

# what about genes unique to Treated_vs_Initial and Untreated_vs_Initial?

unique_sig_TreatedvInitial <- Treated_vs_Initial.p %>% filter(abs(lpv_Treated_vs_Initial) >= 1) %>% anti_join(., Untreated_vs_Initial.p) 
str(unique_sig_TreatedvInitial) #0 genes

unique_sig_Untreated_vs_Initial <- Untreated_vs_Initial.p %>% filter(abs(lpv_Untreated_vs_Initial) >= 1) %>% anti_join(., Treated_vs_Initial.p) 
str(unique_sig_Untreated_vs_Initial )#4 genes


# next look at Treated_vs_Untreated and Untreated_vs_Initial

Treated_vs_Untreated.p %>%
  inner_join(Untreated_vs_Initial.p, by = "gene") %>%
  filter(abs(lpv_Treated_vs_Untreated) >= 1 & abs(lpv_Untreated_vs_Initial) >= 1) %>%
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> TreatedvUntreated_vs_UntreatedvInitial
str(TreatedvUntreated_vs_UntreatedvInitial)
#13 shared genes

# what about genes unique to Treated_vs_Untreated and Untreated_vs_Initial?

unique_sig_TreatedvUntreated <- Treated_vs_Untreated.p %>% filter(abs(lpv_Treated_vs_Untreated) >= 1) %>% anti_join(., Untreated_vs_Initial.p) 
str(unique_sig_TreatedvUntreated) #0 genes

unique_sig_UntreatedvInitial <- Untreated_vs_Initial.p %>% filter(abs(lpv_Untreated_vs_Initial) >= 1) %>% anti_join(., Treated_vs_Untreated.p) 
str(unique_sig_UntreatedvInitial )#0 genes

Untreated_vs_Initial.p %>% 
  filter(abs(lpv_Untreated_vs_Initial) >= 1) %>% 
  select(gene) -> Untreated_vs_Initial_sig

Treated_vs_Initial.p %>% 
  filter(abs(lpv_Treated_vs_Initial) >= 1)  %>% 
  select(gene)-> Treated_vs_Initial_sig

Treated_vs_Untreated.p %>% 
  filter(abs(lpv_Treated_vs_Untreated) >= 1)  %>% 
  select(gene)-> Treated_vs_Untreated_sig
  
pairwise=list("Untreated_vs_Initial"=Untreated_vs_Initial_sig, "Treated_vs_Initial"=Treated_vs_Initial_sig,"Treated_vs_Untreated"=Treated_vs_Untreated_sig)

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

# zero



