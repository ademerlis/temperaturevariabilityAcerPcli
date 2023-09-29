#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### DESEQ IMPORT ####
load("RData_files/pvals.RData")

load("RData_files/control0_control29_lpv.RData")
control0_control29.p <- control0_control29.p %>% dplyr::rename("lpv_c0c29" = lpv) #40,443 genes

load("RData_files/variable0_control0_lpv.RData")
variable0_control0.p <- variable0_control0.p %>% dplyr::rename("lpv_v0c0" = lpv) #20,952 genes

load("RData_files/variable0_variable29_lpv.RData")
variable0_variable29.p <- variable0_variable29.p %>% dplyr::rename("lpv_v0v29" = lpv) #36,731 genes

load("RData_files/variable29_control29_lpv.RData")
variable29_control29.p <- variable29_control29.p %>% dplyr::rename("lpv_v29c29" = lpv) #30,233 genes

#### DEG MATCHING ####

# These sections of code do several things: 1) join common DEGs across experiments with -log10(padj), 
#2) filter by 0.1 padj cutoff (log10(0.1)=1), 3)  adds gene annotations, and 4) then pulls on corresponding KOG classes

# NOTE: normally a -log10 transformation of a p-value between 0 and 1 would not result in a negative number. However, when the above files 
#were generated, (i.e. control0_control29.p), the stat of the directionality of the expression pattern (i.e. negative or positive) was encoded
#in the lpv value manually. The exact lines of code for the transformations were:
#variable29_control29.p$lpv=-log(source[,"padj"],10)
#variable29_control29.p$lpv[source$stat<0]=variable29_control29.p$lpv[source$stat<0]*-1
#This is why the absolute value of lpv is taken for obtaining genes of significance >= 1 (which corresponds to alpha = 0.1)

# first see if there are any shared genes with control0vs29 and variable0vs29
control0_control29.p %>%
  inner_join(variable0_variable29.p, by = "gene") %>%
  filter(abs(lpv_c0c29) >= 1 & abs(lpv_v0v29) >= 1) %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> c0c29_v0v29 #there are ~7,527 shared genes (cut-off is p-adj of 0.1)

# what about genes unique to c0c29 and v0v29?

str(control0_control29.p) #40,443 genes
str(variable0_variable29.p) #36,731 genes

str(control0_control29.p %>% filter(abs(lpv_c0c29) >= 1)) #16,022 genes significant at level of padj 0.1 
str(variable0_variable29.p %>% filter(abs(lpv_v0v29) >= 1)) #10,540 genes significant at level of padj 0.1 

unique_sig_C0C29 <- control0_control29.p %>% filter(abs(lpv_c0c29) >= 1) %>% anti_join(., variable0_variable29.p) #273 genes

unique_sig_V0V29 <- variable0_variable29.p %>% filter(abs(lpv_v0v29) >= 1) %>% anti_join(., control0_control29.p) #0 genes

unique_C0C29 <- anti_join(control0_control29.p, variable0_variable29.p) #3712 genes


# next look at variable0vscontrol0 and variable29vscontrol29
variable0_control0.p %>%
  inner_join(variable29_control29.p, by = "gene") %>%
  filter(abs(lpv_v0c0) >= 1 & abs(lpv_v29c29) >= 1) %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "~/OneDrive - University of Miami/NOAA ERL/stress hardening 2022/gene expression/Acervicornis_annotatedTranscriptome/Acervicornis_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> v0c0_v29c29 #47 shared genes


#### KOG MATCHING ####

# filtering and summarizing DEGs by KOG class for high-level comparisons
v0c0_v29c29 %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v0c0 >= 1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v0c0_up" = n) -> KOG_v0c0_up

v0c0_v29c29 %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v0c0 <= -1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v0c0_down" = n) -> KOG_v0c0_down

v0c0_v29c29 %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v29c29 >= 1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v29c29_up" = n) -> KOG_v29c29_up

v0c0_v29c29 %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v29c29 <= -1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v29c29_down" = n) -> KOG_v29c29_down

c0c29_v0v29 %>% 
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_c0c29 >= 1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "c0c29_up" = n) -> KOG_c0c29_up

c0c29_v0v29 %>% 
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_c0c29 <= -1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "c0c29_down" = n) -> KOG_c0c29_down

c0c29_v0v29 %>% 
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v0v29 >= 1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v0v29_up" = n) -> KOG_v0v29_up

c0c29_v0v29 %>% 
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_v0v29 <= -1) %>%
  dplyr::count(KOG) %>%
  dplyr::rename("KOG" = KOG, "v0v29_down" = n) -> KOG_v0v29_down

# 1) v0c0_v29c29 KOG class sums in a single dataframe

#v0c0 vs. v29c29
KOG_v0c0_up %>%
  inner_join(KOG_v0c0_down, by = "KOG") %>%
  inner_join(KOG_v29c29_down, by = "KOG") %>% 
  inner_join(KOG_v29c29_up, by = "KOG") -> KOG_v0c0_v29c29_match

#c0c29 vs v0v29
KOG_c0c29_down %>%
  inner_join(KOG_c0c29_up, by = "KOG") %>%
  inner_join(KOG_v0v29_up, by = "KOG") %>%
  inner_join(KOG_v0v29_down, by = "KOG") -> KOG_c0c29_v0v29_match

# melting dataframe for plotting
KOG_v0c0_v29c29_match %>% 
  melt(id = "KOG") %>% 
  dplyr::rename(comparison = variable, sum = value) -> KOG_v0c0_v29c29_melt

KOG_c0c29_v0v29_match %>% 
  melt(id = "KOG") %>% 
  dplyr::rename(comparison = variable, sum = value) -> KOG_c0c29_v0v29_melt

# creating a custom color palette
colorCount = length(unique(KOG_v0c0_v29c29_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_v0c0_v29c29_sum <- ggplot(KOG_v0c0_v29c29_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
KOG_v0c0_v29c29_sum
#ggsave("common genes KOG v0c0_v29c29 abundance.pdf", plot= KOG_v0c0_v29c29_sum, width=8, height=6, units="in", dpi=300)


# 2) c0c29_v0v29 plotting

# creating a custom color palette
colorCount = length(unique(KOG_c0c29_v0v29_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_c0c29_v0v29_sum <- ggplot(KOG_c0c29_v0v29_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
KOG_c0c29_v0v29_sum
#ggsave("common genes KOG c0c29_v0v29 abundance.pdf", plot= KOG_c0c29_v0v29_sum, width=12, height=6, units="in", dpi=300)


#### VENN DIAGRAMS ####

# first creating a set of up/downregulated DEGs by v0c0 and v29c29
v0c0_v29c29 %>%
  filter(lpv_v0c0 >= 1) %>%
  pull(gene) -> v0c0_up

v0c0_v29c29 %>%
  filter(lpv_v0c0 <= -1) %>%
  pull(gene) -> v0c0_down

v0c0_v29c29 %>%
  filter(lpv_v29c29 >= 1) %>%
  pull(gene) -> v29c29_up

v0c0_v29c29 %>%
  filter(lpv_v29c29 <= -1) %>%
  pull(gene) -> v29c29_down

# then creating a second set for c0c29_v0v29
c0c29_v0v29 %>%
  filter(lpv_c0c29 >= 1) %>%
  pull(gene) -> c0c29_up

c0c29_v0v29 %>%
  filter(lpv_c0c29 <= -1) %>%
  pull(gene) -> c0c29_down

c0c29_v0v29 %>%
  filter(lpv_v0v29 >= 1) %>%
  pull(gene) -> v0v29_up

c0c29_v0v29 %>%
  filter(lpv_v0v29 <= -1) %>%
  pull(gene) -> v0v29_down


# v0c0_v29c29
venn_v0c0v29c29=venn.diagram(
  x = list("V0/C0 up"=v0c0_up, "V0/C0 down"=v0c0_down,"V29/C29 up"=v29c29_up, "V29/C29 down"=v29c29_down),
  category.names = c("V0/C0 up", "V0/C0 down", "V29/C29 up", "V29/C29 down"),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="venn_v0c0v29c29.pdf", height=10, width=12)
grid.draw(venn_v0c0v29c29)
dev.off()

# c0c29_v0v29
venn_c0c29_v0v29=venn.diagram(
  x = list("C0/C29 up"=c0c29_up, "C0/C29 down"=c0c29_down,"V0/V29 up"=v0v29_up, "V0/V29 down"=v0v29_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="venn_c0c29_v0v29.pdf", height=10, width=12)
grid.draw(venn_c0c29_v0v29)
dev.off()
