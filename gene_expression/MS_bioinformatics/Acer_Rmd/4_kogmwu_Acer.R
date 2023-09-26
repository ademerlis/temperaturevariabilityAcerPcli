#### packages ####

#install.packages("KOGMWU")
library(KOGMWU)


#### pairwise treatments (lpv) ####

# loading KOG annotations
gene2kog=read.table("~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2kogClass.tab",sep="\t", fill=T) #iso2kogClass.tab not iso2kogClass1.tab because that file has an "error" when you try to view it using the terminal
head(gene2kog)

#setwd("OneDrive - University of Miami/GitHub/Ch2_temperaturevariability2023/gene_expression/MS_bioinformatics/Acer_Rmd")
load("RData_files/pvals.RData")

control0_control29=load('RData_files/control0_control29_lpv.RData')
control0_control29 # names of datasets in the package
lpv.control0_control29=kog.mwu(control0_control29.p,gene2kog) 
lpv.control0_control29 

variable0_control0=load('RData_files/variable0_control0_lpv.RData')
variable0_control0 # names of datasets in the package
lpv.variable0_control0=kog.mwu(variable0_control0.p,gene2kog) 
lpv.variable0_control0

variable29_control29=load('RData_files/variable29_control29_lpv.RData')
variable29_control29 # names of datasets in the package
lpv.variable29_control29=kog.mwu(variable29_control29.p,gene2kog) 
lpv.variable29_control29

variable0_variable29=load('RData_files/variable0_variable29_lpv.RData')
variable0_variable29 # names of datasets in the package
lpv.variable0_variable29=kog.mwu(variable0_variable29.p,gene2kog) 
lpv.variable0_variable29


# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("control0_control29"=lpv.control0_control29,"variable0_control0"=lpv.variable0_control0,"variable0_variable29"=lpv.variable0_variable29,"variable29_control29"=lpv.variable29_control29))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_Acer_host_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval) 

# creating a pub-ready corr plot
pdf(file="KOG_intervention_Acer_corr_lpv.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="control0_control29",y="variable0_control0",ktable)
corrPlot(x="control0_control29",y="variable0_variable29",ktable)
corrPlot(x="control0_control29",y="variable29_control29",ktable)

corrPlot(x="variable0_control0",y="variable0_variable29",ktable)
corrPlot(x="variable0_control0",y="variable29_control29",ktable)

corrPlot(x="variable0_variable29",y="variable29_control29",ktable)

dev.off()


#### pairwise treatments (fc) ####

variable0_control0=load('RData_files/variable0_control0_fc.RData')
variable0_control0 # names of datasets in the package
fc.variable0_control0=kog.mwu(variable0_control0.fc,gene2kog) 
fc.variable0_control0 

control0_control29=load('RData_files/control0_control29_fc.RData')
control0_control29 # names of datasets in the package
fc.control0_control29=kog.mwu(control0_control29.fc,gene2kog) 
fc.control0_control29

variable0_variable29=load('RData_files/variable0_variable29_fc.RData')
variable0_variable29 # names of datasets in the package
fc.variable0_variable29=kog.mwu(variable0_variable29.fc,gene2kog) 
fc.variable0_variable29

variable29_control29=load('RData_files/variable29_control29_fc.RData')
variable29_control29 # names of datasets in the package
fc.variable29_control29=kog.mwu(variable29_control29.fc,gene2kog) 
fc.variable29_control29

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("variable0_control0"=fc.variable0_control0,"control0_control29"=fc.control0_control29,
                                "variable0_variable29"=fc.variable0_variable29,"variable29_control29"=fc.variable29_control29))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_Acer_host_fc.pdf")
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()
#this wasn't working so I manually saved it as a PDF by right clicking it in the R plot panel

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_Acer_host_corr_fc.pdf", width=10, height=10)
par(mfrow=c(2,3))
corrPlot(x="control0_control29",y="variable0_control0",ktable)
corrPlot(x="control0_control29",y="variable0_variable29",ktable)
corrPlot(x="control0_control29",y="variable29_control29",ktable)

corrPlot(x="variable0_control0",y="variable0_variable29",ktable)
corrPlot(x="variable0_control0",y="variable29_control29",ktable)

corrPlot(x="variable0_variable29",y="variable29_control29",ktable)

dev.off()
