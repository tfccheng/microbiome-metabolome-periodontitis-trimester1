library(mia)
library(miaViz)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(scater)
library(ggsignif)
library(ggpubr)
library(scales)
set.seed(123)
setwd("")

#load prefiltered data
tse <- readRDS("tse.silva.walpha.rds")
colData(tse)$Sample <- factor(colData(tse)$Sample, levels = c("SBP", "MSA", "MST"))
tse <- transformCounts(tse, method = "relabundance")


####by sample####
#PERMANOVA
permanova.bray.bysample <- adonis2(t(assay(tse,"relabundance")) ~ Sample,
                                        by = "margin", 
                                        data = colData(tse),
                                        method = "bray",
                                        permutations = 999)
permanova.bray.bydisease <- adonis2(t(assay(tse,"relabundance")) ~ Disease,
                                   by = "margin", 
                                   data = colData(tse),
                                   method = "bray",
                                   permutations = 999)

# Perform PCoA with Bray Curtis dissimilarity

tse <- runMDS(tse, FUN = vegan::vegdist, 
                        method = "bray", name = "PCoA_BC", exprs_values = "relabundance")
p.pcoa.bray.bysample <- plotReducedDim(tse, "PCoA_BC", 
                              colour_by = "Sample", shape_by = "Disease",
                              point_size =5)
e <- attr(reducedDim(tse, "PCoA_BC"), "eig")
rel_eig <- e/sum(e[e>0])  
p.pcoa.bray.bysample <- p.pcoa.bray.bysample + 
  labs(subtitle = paste0("Bray-Curtis; P: ", permanova.bray.bysample$`Pr(>F)`),
       color = "Sample",
       x = paste("PCoA1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))+
  scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1) 
p.pcoa.bray.bysample


# Perform NMDS with Bray Curtis dissimilarity
tse <- runNMDS(tse, FUN = vegan::vegdist, 
                    method = "bray", name = "NMDS_BC", exprs_values = "relabundance")
p.nmds.bray.bysample <- plotReducedDim(tse, "NMDS_BC", 
                                       colour_by = "Sample", shape_by = "Disease",
                              point_size =5)
p.nmds.bray.bysample <- p.nmds.bray.bysample + 
  labs(x = "NMDS1", y = "NMDS2")+
  labs(subtitle = paste0("Bray-Curtis; P: ", permanova.bray.bysample$`Pr(>F)`),
       color = "Sample")+
  scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1)
p.nmds.bray.bysample
attr(reducedDim(tse, "NMDS_BC"), "Stress") #NMDS stress value

ggpubr::ggarrange(p.pcoa.bray.bysample, p.nmds.bray.bysample,
                  nrow = 1, ncol = 2, 
                  common.legend = TRUE, legend = "right")

##PERMANOVA subsample by disease Bray Curtis
tse.subsamples <- list()
permanova.bray.subsample.list <- list()
i <- 1

for (n in c("SBP", "MSA", "MST")){
  tse.subsample <- tse[ , tse$Sample %in% n]
  permanova.bray.subsample <- adonis2(t(assay(tse.subsample,"relabundance")) ~ Disease,
                                 by = "margin", 
                                 data = colData(tse.subsample),
                                 method = "bray",
                                 permutations = 999)
  permanova.bray.subsample.list[[n]] <- permanova.bray.subsample}



