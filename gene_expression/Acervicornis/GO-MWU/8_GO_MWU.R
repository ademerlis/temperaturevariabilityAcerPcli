
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to 
#identify GO categories that are significantly enriches with either up- or down-regulated genes. 
#The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical 
#"GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure.

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or
#down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. 
#Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it;
#"good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). 
#For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################

source("gomwu.functions.R")
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml 
goAnnotations="Acropora_iso2go.tab" # two-column, tab-delimited, one line per gene, 
#multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.


#### fold change files ####

# Edit these to match your data file names: 

input="Treated_vs_Initial_fc.csv" 
# BP 802 
# MF 117
# CC 156

input="Untreated_vs_Initial_fc.csv" 
# BP 651
# MF 105
# CC 149

input="Treated_vs_Untreated_fc.csv" 
#BP 160
# MF 15
# CC 89

goDivision="BP" # either MF, or BP, or CC
goDivision="MF" # either MF, or BP, or CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, 
           #kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# is there a way to not have to re-run gomwuStats but just use the files that have already been created to validate the number of sig GO terms?

library(tidyverse)

MWU_result_BP_TreatedvInitial=read.csv("MWU_BP_Treated_vs_Initial_fc.csv", sep = "")
MWU_result_BP_TreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #802 terms = that matches above yay

MWU_result_BP_TreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #543 terms

MWU_result_BP_TreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  ggplot(., aes(x=delta.rank,y=nseqs)) +
  geom_point() # highest and lowest rankings have the smallest number of genes. I guess that makes sense because the log-2-fold change of most genes is going to be closer to zero? 

allgenes_BP_TreatedvInitial=read.csv("BP_Treated_vs_Initial_fc.csv", sep="")
MWU_result_BP_TreatedvInitial %>% 
  filter(p.adj < 0.05) %>%  #543 rows
  left_join(., allgenes_BP_TreatedvInitial, by = c("name", "term")) %>% view()

MWU_result_MF_TreatedvInitial=read.csv("MWU_MF_Treated_vs_Initial_fc.csv", sep = "")

MWU_result_MF_TreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #117 terms = that matches above yay

MWU_result_MF_TreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #94 terms

MWU_result_CC_TreatedvInitial=read.csv("MWU_CC_Treated_vs_Initial_fc.csv", sep = "")

MWU_result_CC_TreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #156 terms = that matches above yay

MWU_result_CC_TreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #111 terms

MWU_result_BP_TreatedvUntreated=read.csv("MWU_BP_Treated_vs_Untreated_fc.csv", sep = "")

MWU_result_BP_TreatedvUntreated %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #160 terms = that matches above yay

MWU_result_BP_TreatedvUntreated %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #114 terms

MWU_result_MF_TreatedvUntreated=read.csv("MWU_MF_Treated_vs_Untreated_fc.csv", sep = "")

MWU_result_MF_TreatedvUntreated %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #15 terms = that matches above yay

MWU_result_MF_TreatedvUntreated %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #10 terms

MWU_result_CC_TreatedvUntreated=read.csv("MWU_CC_Treated_vs_Untreated_fc.csv", sep = "")

MWU_result_CC_TreatedvUntreated %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #89 terms = that matches above yay

MWU_result_CC_TreatedvUntreated %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #74 terms

MWU_result_BP_UntreatedvInitial=read.csv("MWU_BP_Untreated_vs_Initial_fc.csv", sep = "")

MWU_result_BP_UntreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #651 terms = that matches above yay

MWU_result_BP_UntreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #432 terms

MWU_result_MF_UtreatedvInitial=read.csv("MWU_MF_Untreated_vs_Initial_fc.csv", sep = "")

MWU_result_MF_UtreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #105 terms = that matches above yay

MWU_result_MF_UtreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #70 terms

MWU_result_CC_UtreatedvInitial=read.csv("MWU_CC_Untreated_vs_Initial_fc.csv", sep = "")

MWU_result_CC_UtreatedvInitial %>% 
  filter(p.adj < 0.1) %>% 
  summarise(count = n()) #149 terms = that matches above yay

MWU_result_CC_UtreatedvInitial %>% 
  filter(p.adj < 0.05) %>% 
  summarise(count = n()) #104 terms


# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". 
                  #Specify -log(0.05,10) for log p-value (lpv) datasets, and 1 for fold change (fc) datasets. 
                  #Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module 
                  #(all non-zero genes = "good genes").
                  absValue=1,
                  #absValue=0.001,
                  # level1=1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  # level1=0.1,
                  level1=0.05,
                  # level2=0.1, # FDR cutoff to print in regular (not italic) font.
                  # level2=0.05,
                  level2=0.01,
                  # level3=0.05, # FDR cutoff to print in large bold font.
                  # level3=0.01,
                  level3=0.001,
                  txtsize=1.2,    # decrease to fit more on one page, or increase 
                  #(after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results
write.csv(results, file = "TreatedvsUntreated_CC_fc_05.csv")


#### Fisher exact test ####

# I want to just show GO categories over-represented among the genes that have 1 as their measure.
# First I can take the list of DGEs for each contrast and make the abs(lpv) > 1.3 (which corresponds to significance level of 0.05) be equal to 1, and then all other genes = 0.

library(tidyverse)

treated_vs_initial=read.csv("Treated_vs_Initial_lpv.csv")

treated_vs_initial %>% 
  mutate(lpv = case_when(abs(lpv) > 1.3 ~ 1,
                         abs(lpv) <= 1.3 ~ 0)) -> treated_vs_initial

treated_vs_initial %>% 
  group_by(lpv) %>% 
  summarise(count = n()) #5736 genes = this matches results yay

treated_vs_initial %>% 
  write_csv("Treated_vs_Initial_lpv_fisher.csv")

input="Treated_vs_Initial_lpv_fisher.csv" 

goDivision="BP" # either MF, or BP, or CC
goDivision="MF" # either MF, or BP, or CC
#goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, 
           #kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". 
                  #Specify -log(0.05,10) for log p-value (lpv) datasets, and 1 for fold change (fc) datasets. 
                  #Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module 
                  #(all non-zero genes = "good genes").
                  #absValue=1,
                  absValue=0.001,
                  # level1=1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level1=0.1,
                  #level1=0.05,
                  level2=0.1, # FDR cutoff to print in regular (not italic) font.
                  # level2=0.05,
                  #level2=0.01,
                  level3=0.05, # FDR cutoff to print in large bold font.
                  # level3=0.01,
                  #level3=0.001,
                  txtsize=1.2,    # decrease to fit more on one page, or increase 
                  #(after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5 # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

#### log-padj value files ####

# Edit these to match your data file names: 

# two columns of comma-separated values: gene id, continuous measure of significance. 

input="Treated_vs_Initial_lpv.csv" 
# MF 151 GO terms at 10% FDR 
# BP 869
# CC 185

input="Untreated_vs_Initial_lpv.csv" 
# MF 129
# BP 684
# CC 176

input="Treated_vs_Untreated_lpv.csv" 
# MF 24
# BP 218 
# CC 100

goDivision="MF" # either MF, or BP, or CC
goDivision="BP" # either MF, or BP, or CC
goDivision="CC" # either MF, or BP, or CC

# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, 
#go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, 
           #kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". 
                  #Specify -log(0.05,10) for log p-value (lpv) datasets, and 1 for fold change (fc) datasets. 
                  #Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module 
                  #(all non-zero genes = "good genes").
                  # absValue=1,
                  #absValue=0.001,
                  # level1=1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  # level1=0.1,
                  level1=0.05,
                  # level2=0.1, # FDR cutoff to print in regular (not italic) font.
                  # level2=0.05,
                  level2=0.01,
                  # level3=0.05, # FDR cutoff to print in large bold font.
                  # level3=0.01,
                  level3=0.001,
                  txtsize=1.2,    # decrease to fit more on one page, or increase 
                  #(after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results
write.csv(results, file = "TreatedvsUntreated_lpv_CC_05.csv")



