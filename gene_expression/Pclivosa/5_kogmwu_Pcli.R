#### packages ####

#install.packages("KOGMWU")
library(KOGMWU)
library(tidyverse)

#### loading KOG annotations ####
gene2kog=read.table("bioinformatics/Pclivosa_iso2kogClass.tab",sep="\t", fill=T) #iso2kogClass.tab not iso2kogClass1.tab because that file has an "error" when you try to view it using the terminal
head(gene2kog)
gene2kog %>% 
  rename(gene = V1, KOG = V2) -> gene2kog_table

read.table(file = "bioinformatics/Pclivosa_iso2geneName.tab",
           sep = "\t",
           quote="", fill=FALSE) %>%
  rename(gene = V1,
         annot = V2) -> iso2geneName

#### pairwise treatments (fc) ####

load("RData_files/Treated_vs_Initial_fc.RData") 
fc.Treated_vs_Initial=kog.mwu(Treated_vs_Initial.fc,gene2kog)
#write_csv(fc.Treated_vs_Initial, "KOG_pvalues_fc_TreatedvInitial.csv")

full_join(Treated_vs_Initial.fc, gene2kog_table, by="gene") %>% 
  drop_na() %>% 
  filter(!KOG == "" & !KOG == "Function Unknown") %>% 
  rename(term = KOG) %>% 
  full_join(., fc.Treated_vs_Initial, by = "term") %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  select(gene, annot, lfc, term, nseqs:padj) %>% 
  rename(KOG = term) %>% 
  write_csv("KOGterms_allgenes_pvalues_fc_TreatedvInitial.csv")

fc.Treated_vs_Initial %>% 
  filter(!term == "" & !term == "Function Unknown") -> fc.Treated_vs_Initial

fc.Treated_vs_Initial %>% 
  select(term, padj) %>% 
  rename(padj_treatedVinitial = padj) -> df1

load("RData_files/Treated_vs_Untreated_fc.RData") 
fc.Treated_vs_Untreated=kog.mwu(Treated_vs_Untreated.fc,gene2kog)
#write_csv(fc.Treated_vs_Untreated, "KOG_pvalues_fc_TreatedvUntreated.csv")

full_join(Treated_vs_Untreated.fc, gene2kog_table, by="gene") %>% 
  drop_na() %>% 
  filter(!KOG == "" & !KOG == "Function Unknown") %>% 
  rename(term = KOG) %>% 
  full_join(., fc.Treated_vs_Initial, by = "term") %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  select(gene, annot, lfc, term, nseqs:padj) %>% 
  rename(KOG = term) %>% 
  write_csv("KOGterms_allgenes_pvalues_fc_TreatedvUntreated.csv")

fc.Treated_vs_Untreated %>% 
  filter(!term == "" & !term == "Function Unknown") -> fc.Treated_vs_Untreated

fc.Treated_vs_Untreated %>% 
  select(term, padj) %>% 
  rename(padj_treatedVuntreated = padj) -> df2

load("RData_files/Untreated_vs_Initial.fc.RData") 
fc.Untreated_vs_Initial=kog.mwu(Untreated_vs_Initial.fc,gene2kog)
#write_csv(fc.Untreated_vs_Initial, "KOG_pvalues_fv_UntreatedvInitial.csv")

full_join(Untreated_vs_Initial.fc, gene2kog_table, by="gene") %>% 
  drop_na() %>% 
  filter(!KOG == "" & !KOG == "Function Unknown") %>% 
  rename(term = KOG) %>% 
  full_join(., fc.Treated_vs_Initial, by = "term") %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  select(gene, annot, lfc, term, nseqs:padj) %>% 
  rename(KOG = term) %>%  
  write_csv("KOGterms_allgenes_pvalues_fc_UntreatedvInitial.csv")

fc.Untreated_vs_Initial %>% 
  filter(!term == "" & !term == "Function Unknown") -> fc.Untreated_vs_Initial

fc.Untreated_vs_Initial %>% 
  select(term, padj) %>% 
  rename(padj_untreatedVinitial = padj) -> df3

full_join(df1, df2) %>% 
  full_join(., df3) -> df_padj_fc

df_padj_fc %>%
  mutate(
    Treated_vs_Initial = ifelse(padj_treatedVinitial < 0.05, 1, 0),
    Treated_vs_Untreated = ifelse(padj_treatedVuntreated < 0.05, 1, 0),
    Untreated_vs_Initial = ifelse(padj_untreatedVinitial < 0.05, 1, 0)
  ) %>% 
  mutate(Treated_vs_Initial = as.numeric(Treated_vs_Initial),
         Treated_vs_Untreated = as.numeric(Treated_vs_Untreated),
         Untreated_vs_Initial = as.numeric(Untreated_vs_Initial)) %>% 
  column_to_rownames(var="term") %>% 
  select(Treated_vs_Initial:Untreated_vs_Initial) %>% 
  as.matrix() -> df_transformed_fc


#each fc table has 23 KOG terms. I think this is a setting within the kog.mwu function

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Treated_vs_Initial"=fc.Treated_vs_Initial, "Treated_vs_Untreated"=fc.Treated_vs_Untreated, "Untreated_vs_Initial"=fc.Untreated_vs_Initial))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_Pcli_host_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 
while (!is.null(dev.list()))  dev.off()
#needed to manually save this as a PDF

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval) 


#### pairwise treatments (lpv) ####

load("RData_files/pvals.RData")

load("RData_files/Treated_vs_Initial_lpv.RData") #Treated_vs_Initial.p dataset
lpv.Treated_vs_Initial=kog.mwu(Treated_vs_Initial.p,gene2kog)
#write_csv(lpv.Treated_vs_Initial, "KOG_pvalues_lpv_TreatedVInitial.csv")
lpv.Treated_vs_Initial %>% 
  filter(!term == "" & !term == "Function Unknown") -> lpv.Treated_vs_Initial

lpv.Treated_vs_Initial %>% 
  select(term, padj) %>% 
  rename(padj_treatedVinitial = padj) -> df1

load("RData_files/Treated_vs_Untreated_lpv.RData") #Treated_vs_Untreated.p dataset
lpv.Treated_vs_Untreated=kog.mwu(Treated_vs_Untreated.p,gene2kog)
# write_csv(lpv.Treated_vs_Untreated, "KOG_pvalues_lpv_TreatedvUntreated.csv")
lpv.Treated_vs_Untreated %>% 
  filter(!term == "" & !term == "Function Unknown") -> lpv.Treated_vs_Untreated

lpv.Treated_vs_Untreated %>% 
  select(term, padj) %>% 
  rename(padj_treatedVuntreated = padj) -> df2

load("RData_files/Untreated_vs_Initial_lpv.RData") #Treated_vs_Untreated.p dataset
lpv.Untreated_vs_Initial=kog.mwu(Untreated_vs_Initial.p,gene2kog)
#write_csv(lpv.Untreated_vs_Initial, "KOG_pvalues_lpv_UntreatedvInitial.csv")
lpv.Untreated_vs_Initial %>% 
  filter(!term == "" & !term == "Function Unknown") -> lpv.Untreated_vs_Initial

lpv.Untreated_vs_Initial %>% 
  select(term, padj) %>% 
  rename(padj_untreatedVinitial = padj) -> df3

full_join(df1, df2) %>% 
  full_join(., df3) -> df_padj

df_padj %>%
  mutate(
    Treated_vs_Initial = ifelse(padj_treatedVinitial < 0.05, 1, 0),
    Treated_vs_Untreated = ifelse(padj_treatedVuntreated < 0.05, 1, 0),
    Untreated_vs_Initial = ifelse(padj_untreatedVinitial < 0.05, 1, 0)
  ) %>% 
  mutate(Treated_vs_Initial = as.numeric(Treated_vs_Initial),
         Treated_vs_Untreated = as.numeric(Treated_vs_Untreated),
         Untreated_vs_Initial = as.numeric(Untreated_vs_Initial)) %>% 
  column_to_rownames(var="term") %>% 
  select(Treated_vs_Initial:Untreated_vs_Initial) %>% 
  as.matrix() -> df_transformed

ktable_matrix <- as.matrix(ktable)

#each lpv table has 23 KOG terms. I think this is a setting within the kog.mwu function

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Treated_vs_Initial"=lpv.Treated_vs_Initial, "Treated_vs_Untreated"=lpv.Treated_vs_Untreated, "Untreated_vs_Initial"=lpv.Untreated_vs_Initial))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_Pcli_host_lpv.pdf", width=7, height=8)
pheatmap::pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 

#note: neither of these worked 
# pheatmap::pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, 
#                    cell_fun = function(j, i, x, y, width, height, fill) {
#                      if (abs(ktable_matrix[i, j]) > 500)  {  # Check if the value is significant
#                        rect(x, y, x + width, y + height, col = "black", border = "black")
#                      }
#                    }) 
# 
# pheatmap(as.matrix(ktable), clustering_distance_cols="correlation",
#          cell_fun = function(j, i, x, y, width, height, fill) {
#            if (abs(ktable_matrix[i, j]) > 500)  { 
#            rect(x, y, x + width, y + height, col = "black", border = "red", lwd = 2)
#          }})

while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval) 
