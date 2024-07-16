library(DESeq2)
library(edgeR)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(reshape2)
library(mia)
library(microbiomeMarker)
set.seed(123)
setwd("")
tse <- readRDS("tse.silva.wclinical.rds")
df.clinical.num.impute <- readRDS("D:/Research Data/SZ project/Metabolomics/T1 data/New analysis/Data prepare/df.clinical.num.impute.rds")
colData(tse)$`BMI_T1.x` <- rep(df.clinical.num.impute$BMI, 3)
names(rowData(tse))[1] <- "Kingdom"

####subset by Samples at ASV levels####
tse.subsamples <- list()
pseq.subsamples <- list()
##Data preparation with subsample
for (n in c("SBP", "MSA", "MST")){
  tse.subsample <- tse[ , tse$Sample %in% n]
  tse.subsample <- subsetByPrevalentTaxa(tse.subsample, detection = 1E-5, prevalence = 0.05) 
  pseq.subsample <- makePhyloseqFromTreeSummarizedExperiment(tse.subsample)
  tse.subsample <- transformAssay(tse.subsample, assay.type = "counts", method = "log10", MARGIN = "samples", pseudocount = 0.1) 
  tse.subsamples[[n]] <- tse.subsample
  pseq.subsamples[[n]] <- pseq.subsample
}
##DESeq2
sigtab_DESq2_subsamples <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", 
                      taxa_rank = "none", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples[[1]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_SBP.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples[[1]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_SBP.rds")
write.csv(sigtab_DESq2_subsamples[[2]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MSA.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples[[2]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MSA.rds")
write.csv(sigtab_DESq2_subsamples[[3]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MST.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples[[3]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MST.rds")


sigtab_DESq2_subsamples.BMI <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", confounders = "BMI_T1.x",
                      taxa_rank = "none", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.BMI[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.BMI[[1]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_SBP.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.BMI[[1]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_SBP.BMI.rds")
write.csv(sigtab_DESq2_subsamples.BMI[[2]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MSA.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.BMI[[2]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MSA.BMI.rds")
write.csv(sigtab_DESq2_subsamples.BMI[[3]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MST.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.BMI[[3]], "Figures/DESeq2_FD0.05/sig.ASV_DESeq2_MST.BMI.rds")



####subset by Samples at Species levels####
tse.subsamples.sp <- list()
pseq.subsamples.sp <- list()
##Data preparation with subsample
for (n in c("SBP", "MSA", "MST")){
  tse.subsample <- tse[ , tse$Sample %in% n]
  tse.subsample <- subsetByPrevalentTaxa(tse.subsample, detection = 1E-5, prevalence = 0.05)
  tse.subsample <- agglomerateByRank(tse.subsample, rank = "Species", onRankOnly = T)
  pseq.subsample <- makePhyloseqFromTreeSummarizedExperiment(tse.subsample)
  tse.subsample <- transformAssay(tse.subsample, assay.type = "counts", method = "log10", MARGIN = "samples", pseudocount = 0.1)
  tse.subsamples.sp[[n]] <- tse.subsample
  pseq.subsamples.sp[[n]] <- pseq.subsample
}

##DESeq2
sigtab_DESq2_subsamples.sp <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples.sp[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", 
                      taxa_rank = "Species", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.sp[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.sp[[1]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_SBP.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp[[1]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_SBP.rds")
write.csv(sigtab_DESq2_subsamples.sp[[2]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MSA.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp[[2]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MSA.rds")
write.csv(sigtab_DESq2_subsamples.sp[[3]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MST.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp[[3]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MST.rds")


sigtab_DESq2_subsamples.sp.BMI <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples.sp[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", confounders = "BMI_T1.x",
                      taxa_rank = "Species", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.sp.BMI[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.sp.BMI[[1]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_SBP.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp.BMI[[1]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_SBP.BMI.rds")
write.csv(sigtab_DESq2_subsamples.sp.BMI[[2]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MSA.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp.BMI[[2]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MSA.BMI.rds")
write.csv(sigtab_DESq2_subsamples.sp.BMI[[3]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MST.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.sp.BMI[[3]], "Figures/DESeq2_FD0.05/sig.sp_DESeq2_MST.BMI.rds")




#### Genus levels ####
#subset by Samples at Genus levels
tse.subsamples.genus <- list()
pseq.subsamples.genus <- list()
##Data preparation with subsample
for (n in c("SBP", "MSA", "MST")){
  tse.subsample <- tse[ , tse$Sample %in% n]
  tse.subsample <- subsetByPrevalentTaxa(tse.subsample, detection = 1E-5, prevalence = 0.05)
  tse.subsample <- agglomerateByRank(tse.subsample, rank = "Genus", onRankOnly = T)
  pseq.subsample <- makePhyloseqFromTreeSummarizedExperiment(tse.subsample)
  tse.subsample <- transformAssay(tse.subsample, assay.type = "counts", method = "log10", MARGIN = "samples", pseudocount = 0.1)  
  tse.subsamples.genus[[n]] <- tse.subsample
  pseq.subsamples.genus[[n]] <- pseq.subsample
}
##DESeq2
sigtab_DESq2_subsamples.genus <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples.genus[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", 
                      taxa_rank = "Genus", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.genus[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.genus[[1]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus[[1]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.rds")
write.csv(sigtab_DESq2_subsamples.genus[[2]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus[[2]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.rds")
write.csv(sigtab_DESq2_subsamples.genus[[3]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus[[3]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.rds")


sigtab_DESq2_subsamples.genus.BMI <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples.genus[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", confounders = "BMI_T1.x",
                      taxa_rank = "Genus", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.genus.BMI[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.genus.BMI[[1]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus.BMI[[1]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.BMI.rds")
write.csv(sigtab_DESq2_subsamples.genus.BMI[[2]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus.BMI[[2]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.BMI.rds")
write.csv(sigtab_DESq2_subsamples.genus.BMI[[3]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.genus.BMI[[3]], "Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.BMI.rds")

p_dds_subsamples.genus <- list()
i = 1
for (n in c("SBP", "MSA", "MST")){
  sigtab <- sigtab_DESq2_subsamples.genus[[n]]
  #draw barplot
  p_dds<-ggplot(sigtab, aes(x=reorder(feature, +ef_logFC), y=ef_logFC, 
                            fill=Family)) + 
    theme_bw()+
    geom_bar(stat="identity", color="black", width=0.5, size = 0.25) + coord_flip() +
    labs(x="Genus", y = "log2FC", subtitle = sample.names[i]) +
    geom_hline(yintercept=c(-1,1), linetype="dotted") +
    scale_x_discrete(position = "top") +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.25),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face="italic", color = "black", vjust=0.5),
          axis.line.x = element_line(colour = "black", linewidth = 0.5),
          axis.line.y = element_line(colour = "black", linewidth = 1)) 
  p_dds_subsamples.genus[[n]] <- p_dds
  i = i+1
}
heights_dds.genus <- c()
for (n in c("SBP", "MSA", "MST")){
  heights_dds.genus <- append(heights_dds.genus, 
                        dim(sigtab_DESq2_subsamples.genus[[n]])[1])
}

egg::ggarrange(plots = p_dds_subsamples.genus, ncol = 1, heights = heights_dds.genus)


p_dds_subsamples.genus.BMI <- list()
i = 1
for (n in c("SBP", "MSA", "MST")){
  sigtab <- sigtab_DESq2_subsamples.genus.BMI[[n]]
  #draw barplot
  p_dds<-ggplot(sigtab, aes(x=reorder(feature, +ef_logFC), y=ef_logFC, 
                            fill=Family)) + 
    theme_bw()+
    geom_bar(stat="identity", color="black", width=0.5, size = 0.25) + coord_flip() +
    labs(x="Genus", y = "log2FC", subtitle = sample.names[i]) +
    geom_hline(yintercept=c(-1,1), linetype="dotted") +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.25),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face="italic", color = "black", vjust=0.5),
          axis.line.x = element_line(colour = "black", linewidth = 0.5),
          axis.line.y = element_line(colour = "black", linewidth = 1)) 
  p_dds_subsamples.genus.BMI[[n]] <- p_dds
  i = i+1
}
heights_dds.genus.BMI <- c()
for (n in c("SBP", "MSA", "MST")){
  heights_dds.genus.BMI <- append(heights_dds.genus.BMI, 
                              dim(sigtab_DESq2_subsamples.genus.BMI[[n]])[1])
}

egg::ggarrange(plots = p_dds_subsamples.genus.BMI, ncol = 1, heights = heights_dds.genus.BMI)

#####boxplot####
plots.log10.sig.subsamples.genus <- list()
for (n in c("SBP", "MSA", "MST")){ 
  subsample.sig.genus <- sigtab_DESq2_subsamples.genus[[n]] %>%
    arrange(ef_logFC)
  tse.subsample.sig.genus <- tse.subsamples.genus[[n]][subsample.sig.genus$feature]
  df.log10.subsample.sig.genus <- merge(colData(tse.subsample.sig.genus)[1:3], 
                                        t(assay(tse.subsample.sig.genus, "log10")),
                                        by = 'row.names', all = TRUE) %>%
    tibble::column_to_rownames("Row.names")
  plots.log10.subsample.sig.genus <- ggboxplot(melt(df.log10.subsample.sig.genus), 
                                               x = "variable", y = "value",
                                               color = "Disease", fill = "Disease", alpha = 0.2, rotate = T,
                                               add = "jitter") +
    labs(x ="", y="log10(Counts)") + 
    scale_fill_manual(values=c("#00468B", "#AD002A")) +
    scale_color_manual(values=c("#00468B", "#AD002A")) +
    theme(legend.position = "none") +
    stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", hide.ns = T) 
  plots.log10.sig.subsamples.genus[[n]] <- plots.log10.subsample.sig.genus
}


egg::ggarrange(plots = plots.log10.sig.subsamples.genus, ncol = 1, heights = heights_dds.genus)




#### Family levels ####
#subset by Samples at Family levels
tse.subsamples.family <- list()
pseq.subsamples.family <- list()
##Data preparation with subsample
for (n in c("SBP", "MSA", "MST")){
  tse.subsample <- tse[ , tse$Sample %in% n]
  tse.subsample <- subsetByPrevalentTaxa(tse.subsample, detection = 1E-5, prevalence = 0.05)
  tse.subsample <- agglomerateByRank(tse.subsample, rank = "Family", onRankOnly = T)
  pseq.subsample <- makePhyloseqFromTreeSummarizedExperiment(tse.subsample)
  tse.subsample <- transformAssay(tse.subsample, assay.type = "counts", method = "log10", MARGIN = "samples", pseudocount = 0.1)  
  tse.subsamples.family[[n]] <- tse.subsample
  pseq.subsamples.family[[n]] <- pseq.subsample
}
##DESeq2
sigtab_DESq2_subsamples.family <- list()
for (n in c("SBP", "MSA", "MST")){
  pseq.subsample <- pseq.subsamples.family[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", 
                      taxa_rank = "Family", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.family[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.family[[1]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_SBP.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.family[[1]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_SBP.rds")
write.csv(sigtab_DESq2_subsamples.family[[2]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MSA.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.family[[2]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MSA.rds")
write.csv(sigtab_DESq2_subsamples.family[[3]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MST.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.family[[3]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MST.rds")


sigtab_DESq2_subsamples.family.BMI <- list()
for (n in c("SBP", "MST")){ #no sig feature found in MSA
  pseq.subsample <- pseq.subsamples.family[[n]]
  tax_tab <- tax_table(pseq.subsample) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  sigtab = run_deseq2(pseq.subsample, group = "Disease", confounders = "BMI_T1.x",
                      taxa_rank = "Family", fitType = "parametric", sfType = "poscounts", 
                      p_adjust = "BH", pvalue_cutoff = 0.2)
  sigtab_tax <- dplyr::left_join(data.frame(sigtab@marker_table), tax_tab, by = "feature") 
  sigtab_DESq2_subsamples.family.BMI[[n]] <- sigtab_tax
}

write.csv(sigtab_DESq2_subsamples.family.BMI[[1]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_SBP.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.family.BMI[[1]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_SBP.BMI.rds")

write.csv(sigtab_DESq2_subsamples.family.BMI[[2]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MST.BMI.csv", row.names = F)
saveRDS(sigtab_DESq2_subsamples.family.BMI[[2]], "Figures/DESeq2_FD0.05/sig.family_DESeq2_MST.BMI.rds")
