library(mia)
library(miaViz)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scater)
library(ggsignif)
library(ggpubr)
library(ggsci)
set.seed(123)
setwd("")

#load prefiltered data
pseq.silva.minimal.filter <- readRDS("pseq.silva.minimal.filter.rds")
pseq.gtdb.minimal.filter <- readRDS("pseq.gtdb.minimal.filter.rds")
tab.alpha <-microbiome::alpha(pseq.silva.minimal.filter, index = "all")
tse.silva <- makeTreeSummarizedExperimentFromPhyloseq(pseq.silva.minimal.filter)
colData(tse.silva) <- cbind(colData(tse.silva), tab.alpha) 
tse.silva <- mia::estimateDiversity(tse.silva, 
                                    abund_values = "counts",
                                    index = "faith", 
                                    name = "faith")
saveRDS(tse.silva, "tse.silva.walpha.rds")
colData(tse.silva)$Source <- factor(colData(tse.silva)$Source,
                                    levels = c("SBP.NPD", "SBP.PD",
                                               "MSA.NPD", "MSA.PD",
                                               "MST.NPD", "MST.PD"))
colData(tse.silva)$Disease <- factor(colData(tse.silva)$Disease,
                                    levels = c("NPD", "PD"))
#by source between NPD and PD 
#Method1 all alpha indices
plots.source <- lapply(colnames(colData(tse.silva))[-(1:3)],
                plotColData,
                object = tse.silva,
                x = "Source",
                colour_by = "Source")
plots.source <- lapply(plots.source,"+", theme(axis.text.x = element_text(angle=45,hjust=1)))
#ylabs = c("Observed", "Chao1", "Shannon","Inverse Simpson", "Fisher", "Faith")
comb <- list(c("SBP.NPD", "SBP.PD"), 
             c("MSA.NPD", "MSA.PD"), 
             c("MST.NPD", "MST.PD"))

for (i in 1:23){
  plots.source[[i]] <- plots.source[[i]] +
    labs(x = NULL) +
    geom_signif(comparisons = comb, 
                test= wilcox.test,
                map_signif_level = F) 
}

ggpubr::ggarrange(plotlist = plots.source, nrow = 4, ncol = 6, 
                  common.legend = TRUE, legend = "right")
#method 2
plots.source2 <- list()
i<-1
for (n in c("observed", "chao1", 
            "diversity_shannon","diversity_inverse_simpson", 
            "diversity_fisher", "faith")){
  plots.source2[[n]] <- ggboxplot(as.data.frame(colData(tse.silva)), 
                                  x = "Source", y = n,
                                  color = "Source", fill = "Source", alpha = 0.2,
                                  notch = TRUE, 
                                  add = "jitter") +
    labs(y=ylabs[i]) + 
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    stat_compare_means(method = "wilcox.test", comparisons = comb, label = "p.signif")+
    theme(axis.title.x = element_blank())+ 
    scale_color_npg() +
    scale_fill_npg()
  i = i+1
}
ggpubr::ggarrange(plotlist = plots.source2, nrow = 2, ncol = 3, 
                  common.legend = TRUE, legend = "right")




#by disease between samples
comb1 <- list(c("SBP", "MSA"),
              c("SBP", "MST"), 
              c("MSA", "MST"))
plots.disease <- list()
i <- 1
for (n in c("observed", "chao1", 
            "diversity_shannon","diversity_inverse_simpson", 
            "diversity_fisher", "faith")){
  plots.disease[[n]] <- ggboxplot(as.data.frame(colData(tse.silva)), 
                                  x = "Sample", y = n,
                                  color = "Sample", fill = "Sample", alpha = 0.2,
                                  notch = TRUE, 
                                  facet.by = "Disease", 
                                  palette = c(SBP="#00a087", MSA="#0f77c1", MST ="#e64b35"),
                                  add = "jitter") +
    labs(y=ylabs[i]) + 
    stat_compare_means(method = "wilcox.test", comparisons = comb1, label = "p.signif") +
    theme(axis.title.x = element_blank())+ 
    theme(strip.background = element_blank())
  i = i+1
}
ggpubr::ggarrange(plotlist = plots.disease, nrow = 2, ncol = 3, 
                  common.legend = TRUE, legend = "right")


