library(dplyr)
library(reshape2)
library(mia)
library(psych)
library(corrplot)
library(ggpubr)
library(verification)
library(MLeval)
library(pROC)
library(patchwork)
library(ComplexHeatmap)
library(mixOmics)
library(readxl)
library(writexl)
setwd("")
set.seed(123)
tse.micro <- readRDS("")
taxatable <- rowData(tse.micro)[,taxonomyRanks(tse.micro)] %>% as.data.frame() %>%
  tibble::rownames_to_column("ASV")
sample_meta <- data.frame(colData(tse.micro))
sample_meta$Disease <- factor(sample_meta$Disease)
df.clinical.num.all <- colData(tse.micro)[1:54,-c(1:5)] %>%
  as.data.frame() %>%
  select_if(is.numeric) %>%
  dplyr::select(-age_S)
num_labels <- c("Teeth No.", "FMPS", "BOP", "PD>=4mm",  "CAL 1-2mm",             
                "Age", "Gestational week", "Menarche age", 
                "Menstrual days", "Menstrual cycle", "Height", "Weight", "Current gestational week", 
                "Pregnancy times", "Parity times", "Abortion times", "BMI", "SBP", "DBP", "WBC", "CRP",        
                "UMA", "UCr", "UACR", "Uric acid", "TBIL", 
                "TBA", "ALT", "AST", "TP", "ALB",         
                "TC", "TG", "HDLc", "LDLc", "Ferritin",         
                "Insulin", "Blood glucose", "HbA1c", "TSH", "FT4", "Thyroglobulin", "Fecal calprotectin")
colnames(df.clinical.num.all) <- num_labels

df.clinical.num <- df.clinical.num.all %>%
  dplyr::select(-Age, -Height, -Weight, -`Pregnancy times`, -`Parity times`, -`Abortion times`,
                -`Menstrual days`, -`Menstrual cycle`, -`Current gestational week`)


##Wilcox test
test <- data.frame(matrix(ncol=3,nrow=ncol(df.clinical.num), 
                          dimnames=list(colnames(df.clinical.num), 
                                        c("Shapiro.p", "Wilcox.p", "t.test.p"))))
for (i in rownames(test)){
  shapiro.p <- shapiro.test(df.clinical.num[,i])
  wilcox.p <- wilcox.test(df.clinical.num[,i]~ sample_meta$Disease[1:54])
  ttest.p <- t.test(df.clinical.num[,i]~ sample_meta$Disease[1:54])
  test[i, "Shapiro.p"] <- shapiro.p$p.value
  test[i, "Wilcox.p"] <- wilcox.p$p.value
  test[i, "t.test.p"] <- ttest.p$p.value
}
rownames(test[test$Shapiro.p > 0.05,]) # check normally distributed variables

df.test <- data.frame(df.clinical.num, Disease = sample_meta$Disease[1:54], check.names = F)
df.mean <- df.test %>%
  add_row(Disease = "All", dplyr::summarise(., across(1:34, ~mean(.x, na.rm =T)))) %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~mean(.x, na.rm =T))) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()
df.sd.all <- df.test %>%
  dplyr::summarise(across(1:34, ~sd(.x, na.rm =T))) %>%
  dplyr::mutate(Disease="All", .before = 1)
df.sd <- df.test %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~sd(.x, na.rm =T))) %>%
  rbind(df.sd.all) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()
df.median <- df.test %>%
  add_row(Disease = "All", dplyr::summarise(., across(1:34, ~median(.x, na.rm =T)))) %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~median(.x, na.rm =T))) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()
df.Q1 <- df.test %>%
  add_row(Disease = "All", dplyr::summarise(., across(1:34, ~quantile(.x, 0.25, na.rm =T)))) %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~quantile(.x, 0.25, na.rm =T))) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()
df.Q3 <- df.test %>%
  add_row(Disease = "All", dplyr::summarise(., across(1:34, ~quantile(.x, 0.75, na.rm =T)))) %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~quantile(.x, 0.75, na.rm =T))) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()
df.iqr <- df.test %>%
  add_row(Disease = "All", dplyr::summarise(., across(1:34, ~IQR(.x, na.rm =T)))) %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(across(1:34, ~IQR(.x, na.rm =T))) %>%
  tibble::column_to_rownames("Disease") %>% t() %>% as.data.frame()

test$mean.sd.All <- paste0(sprintf(df.mean$All, fmt = '%#.2f'), 
                           "\u00B1", sprintf(df.sd$All, fmt = '%#.2f'))
test$mean.sd.NPD <- paste0(sprintf(df.mean$NPD, fmt = '%#.2f'), 
                           "\u00B1", sprintf(df.sd$NPD, fmt = '%#.2f'))
test$mean.sd.PD <- paste0(sprintf(df.mean$PD, fmt = '%#.2f'), 
                          "\u00B1", sprintf(df.sd$PD, fmt = '%#.2f'))
test$median.iqr.All <- paste0(sprintf(df.median$All, fmt = '%#.2f'), 
                              "(", sprintf(df.Q1$All, fmt = '%#.2f'), 
                              "\U2013", sprintf(df.Q3$All, fmt = '%#.2f'), ")")
test$median.iqr.NPD <- paste0(sprintf(df.median$NPD, fmt = '%#.2f'), 
                              "(", sprintf(df.Q1$NPD, fmt = '%#.2f'), 
                              "\U2013", sprintf(df.Q3$NPD, fmt = '%#.2f'), ")")
test$median.iqr.PD <- paste0(sprintf(df.median$PD, fmt = '%#.2f'), 
                             "(", sprintf(df.Q1$PD, fmt = '%#.2f'), 
                             "\U2013", sprintf(df.Q3$PD, fmt = '%#.2f'), ")")

write_xlsx(test %>% tibble::rownames_to_column("Parameter"), "Figures/clinical.test.xlsx")



####sPCA, sPLS clinical####
df.clinical.num.non.dent <- df.clinical.num[-c(1:5)]
nipals.tune = nipals(df.clinical.num.non.dent, ncomp = 10)$eig
barplot(nipals.tune, xlab = 'Principal component', ylab = 'Explained variance')
df.clinical.num.impute <- impute.nipals(df.clinical.num.non.dent, ncomp = 3) 
pca.clinical = pca(df.clinical.num.impute, ncomp = 10, center = TRUE, scale = T) 
plot(pca.clinical)
plotIndiv(pca.clinical, group = sample_meta$Disease[1:54], ind.names = FALSE, # plot the samples projected
          legend = TRUE, ellipse = TRUE, title = 'PCA on clinical items, comp 1 - 2')
p.PC1 <- ggplot(data.frame(pca.clinical$variates$X, Disease = sample_meta$Disease[1:54], check.names = F), 
                aes(x=Disease, y=PC1, color = Disease, fill = Disease)) + 
  geom_boxplot(notch=T, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  scale_fill_manual(values=c("#00468B", "#AD002A"), name = "") +
  scale_color_manual(values=c("#00468B", "#AD002A"), name = "") +
  theme_classic() + labs(x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label.x = 1.5, label = "p.format", show.legend = F) 
p.PC2 <- ggplot(data.frame(pca.clinical$variates$X, Disease = sample_meta$Disease[1:54], check.names = F), 
                aes(x=Disease, y=PC2, color = Disease, fill = Disease)) + 
  geom_boxplot(notch=T, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  scale_fill_manual(values=c("#00468B", "#AD002A"), name = "") +
  scale_color_manual(values=c("#00468B", "#AD002A"), name = "") +
  theme_classic() + labs(x = "") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label.x = 1.5, label = "p.format", show.legend = F) 
p.PC12<-ggarrange(p.PC1,p.PC2, ncol = 1, nrow = 2, common.legend = T, legend = "right")



biplot.pca.clinical <- biplot(pca.clinical, cex = 1, ind.names = F,
                              group = sample_meta$Disease[1:54],  # colour by sample class
                              legend.title = '', col.per.group = c("#00468B", "#AD002A"), 
                              title = 'PCA comp 1 - 2')
pl.biplot.pca.clinical <- biplot.pca.clinical +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1)


lay <- rbind(c(1,1,2),
             c(1,1,2))
#plot clinical biplot and barplot of PC12
grid.arrange(pl.biplot.pca.clinical, p.PC12,
             layout_matrix = lay)

plotLoadings(pca.clinical)
test.keepX <- c(seq(5, 29, 5)) # set the number of variable values to be tested

tune.spca.res <- tune.spca(df.clinical.num.impute, ncomp = 3, # generate the first 3 components
                           nrepeat = 50, # repeat the cross-validation 5 times
                           folds = 10, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res)
tune.spca.res$choice.keepX
splsda.clinic <- splsda(df.clinical.num.impute, sample_meta$Disease[1:54], ncomp = 10)
plotIndiv(splsda.clinic , comp = 1:2, 
          group = sample_meta$Disease[1:54], ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')
biplot(splsda.clinic, cex = 0.7, ind.names = F,
       group = sample_meta$Disease[1:54],  # colour by sample class
       legend.title = '', 
       title = 'sPLSDA comp 1 - 2')
perf.splsda.clinic <- perf(splsda.clinic, validation = "Mfold", 
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          progressBar = T, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.clinic, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.clinic$choice.ncomp
list.keepX <- c(seq(3, 29, 2))
tune.splsda.clinic <- tune.splsda(df.clinical.num.impute, sample_meta$Disease[1:54], ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2)
plot(tune.splsda.clinic, col = color.jet(2)) 
tune.splsda.clinic$choice.ncomp$ncomp
tune.splsda.clinic$choice.keepX 
optimal.ncomp <- 2
optimal.keepX <- tune.splsda.clinic$choice.keepX[1:optimal.ncomp]
final.splsda.clinic <- splsda(df.clinical.num.impute, sample_meta$Disease[1:54], 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

plotIndiv(final.splsda.clinic, comp = 1:2, 
          group = sample_meta$Disease[1:54], ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = 'sPLSDA with confidence ellipses')


perf.splsda.clinic <- perf(final.splsda.clinic, 
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = T)
plot(perf.splsda.clinic$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plotLoadings.final.splsda.clinic <- plotLoadings(final.splsda.clinic,comp = 1, contrib = 'max', method = 'median')
df.pld.final.splsda.clinic <- plotLoadings.final.splsda.clinic$X %>%
  tibble::rownames_to_column("parameter") %>% 
  filter(GroupContrib != "tie")
df.clinical.select <- df.pld.final.splsda.clinic %>% dplyr::filter(abs(importance) > 0.1) %>%
  arrange(desc(GroupContrib))
pld.final.splsda.clinic <- ggplot(data=df.pld.final.splsda.clinic, 
                               aes(x=importance, y = reorder(parameter, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Clinical", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


####numeric parameter correlation####
cor.clinical.num <- corr.test(x = df.clinical.num, y= df.clinical.num,
          use = "pairwise",method="spearman",adjust="BH", 
          alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(cor.clinical.num$r,   
         p.mat = cor.clinical.num$p, tl.col = 'black', tl.srt = 45,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig')

r<- cor.clinical.num$r
p <- cor.clinical.num$p %>% as.data.frame()



Heatmap(r, name = "Spearman's \u03C1",
                   cluster_rows = T, cluster_columns = T,
                   row_names_side = "left", column_names_side = "top",
                   cluster_column_slices =F, column_title=NULL,
                   cell_fun = cell_fun, 
                   col = col_fun(9), column_names_rot = 45, 
                   width = ncol(r)*unit(6, "mm"), 
                   height = nrow(r)*unit(6, "mm")) 
print(p.g.SBP)



cor.clinical.num.sig.index <- as.data.frame(cor.clinical.num$p.adj < 0.05)
diag(cor.clinical.num.sig.index) <- FALSE #to prevent the influence of diagonal
not.all.sig.cols <- colSums(cor.clinical.num.sig.index, na.rm = T)
not.all.sig.rows <- rowSums(cor.clinical.num.sig.index, na.rm = T)
cor.clinical.num.r.sig<-cor.clinical.num$r[not.all.sig.rows>0,not.all.sig.cols>0]
cor.clinical.num.padj.sig<-cor.clinical.num$p.adj[not.all.sig.rows>0,not.all.sig.cols>0]
corrplot(cor.clinical.num.r.sig,   
         p.mat = cor.clinical.num.padj.sig, tl.col = 'black', tl.srt = 45, order = 'hclust',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig')
####prepare data
tse.micro.subsamples <- list()
for (n in c("SBP", "MSA", "MST")){
  tse.micro.subsample <- tse.micro[, tse.micro$Sample %in% n] %>%
    subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05)  #filter with 0.05% detection & 5% prevalence
  for (m in c("Genus", "Species")){
    altExp(tse.micro.subsample, m) <- agglomerateByRank(tse.micro.subsample, rank = m, onRankOnly = T)
  }
  tse.micro.subsample <- tse.micro.subsample %>%
    transformSamples(assay_name = "counts", method = "relabundance") %>%
    transformSamples(assay_name = "relabundance", method = "rclr")
  tse.micro.subsamples[[n]] <- tse.micro.subsample
}



sig.genus_DESeq2_SBP <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.rds") %>%
  dplyr::mutate(test = "Disease", .after = 1)
sig.genus_DESeq2_MSA <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.rds") %>%
  dplyr::mutate(test = "Disease", .after = 1)
sig.genus_DESeq2_MST <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.rds") %>%
  dplyr::mutate(test = "Disease", .after = 1)
sig.genus_DESeq2_SBP.BMI <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_SBP.BMI.rds") %>%
  dplyr::mutate(test = "Disease_BMI", .after = 1)
sig.genus_DESeq2_MSA.BMI <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_MSA.BMI.rds") %>%
  dplyr::mutate(test = "Disease_BMI", .after = 1)
sig.genus_DESeq2_MST.BMI <- readRDS("Figures/DESeq2_FD0.05/sig.genus_DESeq2_MST.BMI.rds") %>%
  dplyr::mutate(test = "Disease_BMI", .after = 1)

df.sig.genus.DESeq2.SBP <- dplyr::full_join(sig.genus_DESeq2_SBP[1:6], 
                                      sig.genus_DESeq2_SBP.BMI[1:6], by = "feature") %>%
  dplyr::mutate(enrich_group = ifelse(is.na(enrich_group.x), enrich_group.y, enrich_group.x), .after = 1) %>%
  arrange(desc(enrich_group))
df.sig.genus.DESeq2.MSA <- dplyr::full_join(sig.genus_DESeq2_MSA[1:6], 
                                      sig.genus_DESeq2_MSA.BMI[1:6], by = "feature") %>%
  dplyr::mutate(enrich_group = ifelse(is.na(enrich_group.x), enrich_group.y, enrich_group.x), .after = 1) %>%
  arrange(desc(enrich_group))
df.sig.genus.DESeq2.MST <- dplyr::full_join(sig.genus_DESeq2_MST[1:6], 
                                      sig.genus_DESeq2_MST.BMI[1:6], by = "feature") %>%
  dplyr::mutate(enrich_group = ifelse(is.na(enrich_group.x), enrich_group.y, enrich_group.x), .after = 1) %>%
  arrange(desc(enrich_group)) %>%
  filter(test.x == "Disease") # because BMI adjustment included much more genera for MST, we keep overlapped genera instead 




df.clinical.num.select.perio <- cbind(df.clinical.num[,2:5],df.clinical.num[,df.clinical.select$parameter])
df.clinical.num.group <- data.frame(colnames(df.clinical.num.select.perio), 
                              c(rep("PD", 12), rep("NPD", 12)))
colnames(df.clinical.num.group) <- c("Item", "Group")
####Genus level####

tse.subsamples.g <- list()
for (n in c("SBP", "MSA", "MST")){
  tse.subsample.g <- altExp(tse.micro.subsamples[[n]], "Genus") %>%
    transformAssay(assay.type = "counts", method = "relabundance") %>%
    transformAssay(assay.type = "counts", method = "log10", pseudocount = 0.1) %>%
    transformAssay(assay.type = "counts", method = "z", MARGIN = "features") %>%
    transformAssay(assay.type  = "relabundance", method = "clr", pseudocount = 1E-6) %>%
    transformAssay(assay.type = "relabundance", method = "rclr")
  tse.subsamples.g[[n]] <- tse.subsample.g
}

####SBP correlation with clinical parameters####
df.SBP.g.select <- t(as.data.frame(assay(tse.subsamples.g[[1]], "rclr")))[,df.sig.genus.DESeq2.SBP$feature]
df.SBP.g.select.rel <- t(as.data.frame(assay(tse.subsamples.g[[1]], "relabundance")))[,df.sig.genus.DESeq2.SBP$feature]
df.corr.SBP.g.num <- corr.test(x = df.SBP.g.select.rel,
                                     y = df.clinical.num.select.perio, 
                                     use = "pairwise",method="spearman",adjust="BH", 
                                     alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(df.corr.SBP.g.num$r,   
         p.mat = df.corr.SBP.g.num$p.adj, tl.col = 'black', tl.srt = 45,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig')
write_xlsx(df.corr.SBP.g.num$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.SBP.g.FD0.05.num.p.adj.xlsx")
cl.all <- rbind(df.sig.genus.DESeq2.SBP,
                df.sig.genus.DESeq2.MSA,
                df.sig.genus.DESeq2.MST) %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)
col1 <- setNames(c("#00468B", "#AD002A"), c("NPD", "PD"))
col2 <- setNames(c(paletteer::paletteer_d("ggsci::category20_d3", 
                                          length(unique(cl.all$Class)))), c(unique(cl.all$Class)))
require(RColorBrewer)
col_fun = colorRampPalette(colors = rev(brewer.pal(9,"RdBu")))
level_group <- factor(df.clinical.num.group$Group, levels = c("PD", "NPD"))
cell_fun = function(j, i, x, y, w, h, fill) {
  if(p[i, j] < 0.001 & !is.na(p[i, j])) {
    grid.text("***", x, y)
  } else if(p[i, j] < 0.01& !is.na(p[i, j])) {
    grid.text("**", x, y)
  }
  else if(p[i, j] < 0.05& !is.na(p[i, j])) {
    grid.text("*", x, y)
  }
}


r<- df.corr.SBP.g.num$r
p <- df.corr.SBP.g.num$p.adj %>% as.data.frame()
cl <- df.sig.genus.DESeq2.SBP %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)

ha1 = HeatmapAnnotation(`Enriched Group` = df.clinical.num.group$Group, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.genus.DESeq2.SBP$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.g.SBP <- Heatmap(r, name = "Spearman's \u03C1",
                   top_annotation = ha1, left_annotation = ha2, 
                   cluster_rows = F, cluster_columns = F,
                   row_names_side = "left", column_names_side = "top",
                   row_names_gp = gpar(fontface = "italic"),
                   row_title = "Subgingival plaque", row_title_side = "right",
                   row_title_rot = 90,
                   column_split = level_group,
                   cluster_column_slices =F, column_title=NULL,
                   cell_fun = cell_fun, 
                   col = col_fun(9), column_names_rot = 45, 
                   width = ncol(r)*unit(6, "mm"), 
                   height = nrow(r)*unit(6, "mm")) 
print(p.g.SBP)
####MSA correlation with clinical parameters####
df.MSA.g.select <- t(as.data.frame(assay(tse.subsamples.g[[2]], "rclr")))[,df.sig.genus.DESeq2.MSA$feature]
df.MSA.g.select.rel <- t(as.data.frame(assay(tse.subsamples.g[[2]], "relabundance")))[,df.sig.genus.DESeq2.MSA$feature]
df.corr.MSA.g.num <- corr.test(x = df.MSA.g.select.rel,
                               y = df.clinical.num.select.perio, 
                               use = "pairwise",method="spearman",adjust="BH", 
                               alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(df.corr.MSA.g.num$r,   
         p.mat = df.corr.MSA.g.num$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')
write_xlsx(df.corr.MSA.g.num$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.MSA.g.FD0.05.num.p.adj.xlsx")
r<- df.corr.MSA.g.num$r
p <- df.corr.MSA.g.num$p.adj %>% as.data.frame()
cl <- df.sig.genus.DESeq2.MSA %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)


ha1 = HeatmapAnnotation(`Enriched Group` = df.clinical.num.group$Group, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.genus.DESeq2.MSA$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"),
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.g.MSA <- Heatmap(r, #name = "Spearman's \u03C1",
                   show_heatmap_legend = FALSE,
                   top_annotation = NULL, left_annotation = ha2, 
                   cluster_rows = F, cluster_columns = F,
                   row_names_side = "left", column_names_side = "top",
                   row_names_gp = gpar(fontface = "italic"),
                   row_title = "Saliva", row_title_side = "right",
                   row_title_rot = 90,
                   column_split = level_group,
                   cluster_column_slices =F, column_title=NULL,
                   cell_fun = cell_fun, 
                   col = col_fun(9), column_names_rot = 45, 
                   width = ncol(r)*unit(6, "mm"), 
                   height = nrow(r)*unit(6, "mm")) 
print(p.g.MSA)
####MST  clr data have signicant correlation with clinical parameters####
df.MST.g.select <- t(as.data.frame(assay(tse.subsamples.g[[3]], "rclr")))[,df.sig.genus.DESeq2.MST$feature]
df.MST.g.select.rel <- t(as.data.frame(assay(tse.subsamples.g[[3]], "relabundance")))[,df.sig.genus.DESeq2.MST$feature]
df.MST.g.select.clr <- t(as.data.frame(assay(tse.subsamples.g[[3]], "clr")))[,df.sig.genus.DESeq2.MST$feature]
df.corr.MST.g.num <- corr.test(x = df.MST.g.select.clr,
                               y = df.clinical.num.select.perio, 
                               use = "pairwise",method="spearman",adjust="BH", 
                               alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(df.corr.MST.g.num$r,   
         p.mat = df.corr.MST.g.num$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')
write_xlsx(df.corr.MST.g.num$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.MST.g.clr.FD0.05.num.p.adj.xlsx")
r<- df.corr.MST.g.num$r
p <- df.corr.MST.g.num$p.adj %>% as.data.frame()
cl <- df.sig.genus.DESeq2.MST %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)


ha1 = HeatmapAnnotation(`Enriched Group` = df.clinical.num.group$Group, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.genus.DESeq2.MST$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"),
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.g.MST <- Heatmap(r, name = "Spearman's \u03C1",
                   top_annotation = ha1, left_annotation = ha2, 
                   cluster_rows = F, cluster_columns = F,
                   row_names_side = "left", column_names_side = "top",
                   row_names_gp = gpar(fontface = "italic"),
                   row_title = "Feces", row_title_side = "right",
                   row_title_rot = 90,
                   column_split = level_group,
                   cluster_column_slices =F, column_title=NULL,
                   cell_fun = cell_fun, 
                   col = col_fun(9), column_names_rot = 45, 
                   width = ncol(r)*unit(6, "mm"), 
                   height = nrow(r)*unit(6, "mm")) 

hlist <- p.g.SBP %v% p.g.MSA 
print(p.g.MST)


####plaque fecal microb fecal metab correlation####
df.venn.feces.metab <- read.csv("D:/Research Data/SZ project/Metabolomics/T1 data/New analysis/venn/df.venn.feces.csv")
df.venn.feces.metab.select <- df.venn.feces.metab %>%
  filter(adj.P.Val < 0.2) %>%
  arrange(desc(importance))
df_metab_feces_X.log.pareto <- readRDS("D:/Research Data/SZ project/Metabolomics/T1 data/New analysis/Data prepare/df_metab_feces_X.log.pareto.rds") %>%
  t() %>% as.data.frame() %>% tibble::rownames_to_column("ID")

df.venn.feces.metab_data <- df.venn.feces.metab.select[, c("ID", "identification_pld")] %>%
  dplyr::left_join(df_metab_feces_X.log.pareto, by = "ID") %>%
  dplyr::mutate(label_name = paste0(ID, ": ", identification_pld)) %>%
  dplyr::select(-identification_pld) %>%
  relocate(label_name, .after=ID) 


df.venn.feces.metab_data$label_name <- stringr::str_replace(df.venn.feces.metab_data$label_name, " NA", " ")
df.venn.feces.metab_data$label_name <- stringr::str_replace(df.venn.feces.metab_data$label_name, "MS2 ", "")
df.venn.feces.metab_data$label_name <- stringr::str_replace(df.venn.feces.metab_data$label_name, " \\(\\)", "")
df.venn.feces.metab_data <- df.venn.feces.metab_data %>%
  tibble::column_to_rownames("label_name") %>% 
  dplyr::select(-ID) %>%
  t() %>% as.data.frame() 
####correlation between fecal microb and fecal metab####
df.corr.feces.microb.metab <- corr.test(x = df.MST.g.select.rel,
                             y = df.venn.feces.metab_data, 
                             use = "pairwise",method="spearman",adjust="BH", 
                             alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
write_xlsx(df.corr.feces.microb.metab$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.feces.microb.FD0.05.metab.p.adj.xlsx")
metab_level_group <- factor(df.venn.feces.metab.select$GroupContrib, levels = c("PD", "NPD"))


corrplot(df.corr.feces.microb.metab$r,   
         p.mat = df.corr.feces.microb.metab$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')

r<- df.corr.feces.microb.metab$r
p <- df.corr.feces.microb.metab$p.adj %>% as.data.frame()
cl <- df.sig.genus.DESeq2.MST %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)


ha1 = HeatmapAnnotation(`Enriched Group` = df.venn.feces.metab.select$GroupContrib, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.genus.DESeq2.MST$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"),
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.feces.microb.metab <- Heatmap(r, name = "Spearman's \u03C1",
                 show_heatmap_legend = TRUE,
                 top_annotation = ha1, left_annotation = ha2, 
                 cluster_rows = F, cluster_columns = F,
                 row_names_side = "left", column_names_side = "top",
                 row_names_gp = gpar(fontface = "italic"),
                 row_title = "Feces", row_title_side = "right",
                 row_title_rot = 90,
                 column_split = metab_level_group,
                 cluster_column_slices =F, column_title=NULL,
                 cell_fun = cell_fun.padj, 
                 col = col_fun(9), column_names_rot = 45, 
                 width = ncol(r)*unit(6, "mm"), 
                 height = nrow(r)*unit(6, "mm")) 
print(p.feces.microb.metab)
#clr
df.corr.feces.microb.clr.metab <- corr.test(x = df.MST.g.select.clr,
                                        y = df.venn.feces.metab_data, 
                                        use = "pairwise",method="spearman",adjust="BH", 
                                        alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
write_xlsx(df.corr.feces.microb.clr.metab$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.feces.microb.FD0.05.clr.metab.p.adj.xlsx")


corrplot(df.corr.feces.microb.clr.metab$r,   
         p.mat = df.corr.feces.microb.clr.metab$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')

r<- df.corr.feces.microb.clr.metab$r
p <- df.corr.feces.microb.clr.metab$p.adj %>% as.data.frame()
cl <- df.sig.genus.DESeq2.MST %>%
  dplyr::rename(Genus = feature) %>%
  dplyr::left_join(taxatable, by = "Genus") %>%
  distinct(Genus, .keep_all =TRUE)


ha1 = HeatmapAnnotation(`Enriched Group` = df.venn.feces.metab.select$GroupContrib, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.genus.DESeq2.MST$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"),
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.feces.microb.clr.metab <- Heatmap(r, name = "Spearman's \u03C1",
                                show_heatmap_legend = TRUE,
                                top_annotation = ha1, left_annotation = ha2, 
                                cluster_rows = F, cluster_columns = F,
                                row_names_side = "left", column_names_side = "top",
                                row_names_gp = gpar(fontface = "italic"),
                                row_title = "Feces", row_title_side = "right",
                                row_title_rot = 90,
                                column_split = metab_level_group,
                                cluster_column_slices =F, column_title=NULL,
                                cell_fun = cell_fun.padj, 
                                col = col_fun(9), column_names_rot = 45, 
                                width = ncol(r)*unit(6, "mm"), 
                                height = nrow(r)*unit(6, "mm")) 

print(p.feces.microb.clr.metab)

####scatter correaltion between Coprococcus and urobilin
options(ggpubr.parse_aes = FALSE) #let ggpubr to handle non-standard column names
df.feces.g.metab <- cbind(df.MST.g.select.rel,
                          df.venn.feces.metab_data)
p_Co_urobilin<-ggscatter(df.feces.g.metab, x = "Coprococcus", 
                y = "neg-M593T177: L-Urobilin (HMDB0004159)", 
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line", 
                add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                xlab = "Coprococcus\n(relative aboundance)", ylab = "neg-M593T177: L-Urobilin\nlog10(Pareto scaling)") 
p_Co_kyna1 <- ggscatter(df.feces.g.metab, x = "Coprococcus", 
                y = "pos-M190T161: Kynurenic acid (HMDB0000715)", 
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line", 
                add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                xlab = "Coprococcus\n(relative aboundance)", ylab = "pos-M190T161: Kynurenic acid\nlog10(Pareto scaling)")
p_Co_kyna2 <- ggscatter(df.feces.g.metab, x = "Coprococcus", 
                y = "neg-M188T159: Kynurenic acid (HMDB0000715)", 
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line", 
                add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                xlab = "Coprococcus\n(relative aboundance)", ylab = "neg-M188T159: Kynurenic acid\nlog10(Pareto scaling)")
p_Co_D_lactic <- ggscatter(df.feces.g.metab, x = "Coprococcus", 
                y = "neg-M89T138: D-Lactic acid (HMDB0001311)", 
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line", 
                add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                xlab = "Coprococcus\n(relative aboundance)", ylab = "neg-M89T138: D-Lactic acid\nlog10(Pareto scaling)")
p_Co_eriod<-ggscatter(df.feces.g.metab, x = "Coprococcus", 
                         y = "pos-M289T194: Eriodictyol (HMDB0005810)", 
                         color = "black", shape = 21, size = 3, # Points color, shape and size
                         add = "reg.line", 
                         add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                         xlab = "Coprococcus\n(relative aboundance)", ylab = "pos-M289T194: Eriodictyol\nlog10(Pareto scaling)") 

p_Lachno_urobilin<-ggscatter(df.feces.g.metab, x = "Lachnoclostridium", 
                         y = "neg-M593T177: L-Urobilin (HMDB0004159)", 
                         color = "black", shape = 21, size = 3, # Points color, shape and size
                         add = "reg.line", 
                         add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                         xlab = "Lachnoclostridium\n(relative aboundance)", ylab = "neg-M593T177: L-Urobilin\nlog10(Pareto scaling)") 

p_Lachno_hypotau<-ggscatter(df.feces.g.metab, x = "Lachnoclostridium", 
                             y = "pos-M110T101: Hypotaurine (HMDB0000965)", 
                             color = "black", shape = 21, size = 3, # Points color, shape and size
                             add = "reg.line", 
                             add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                             conf.int = TRUE, # Add confidence interval
                             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                             cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                             xlab = "Lachnoclostridium\n(relative aboundance)", ylab = "pos-M110T101: Hypotaurine\nlog10(Pareto scaling)") 

p_Lachno_xanth<-ggscatter(df.feces.g.metab, x = "Lachnoclostridium", 
                            y = "pos-M153T101_2: Xanthine (HMDB0000292)", 
                            color = "black", shape = 21, size = 3, # Points color, shape and size
                            add = "reg.line", 
                            add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                            conf.int = TRUE, # Add confidence interval
                            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                            cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                            xlab = "Lachnoclostridium\n(relative aboundance)", ylab = "pos-M153T101_2: Xanthine\nlog10(Pareto scaling)") 

ggarrange(p_Co_urobilin, p_Co_kyna1, p_Co_kyna2, p_Co_D_lactic, p_Co_eriod,
          p_Lachno_urobilin, p_Lachno_hypotau, p_Lachno_xanth, ncol = 4, nrow = 2)



df.corr.feces.microb.metab.sp <- corr.test(x = df.MST.sp.select.rel,
                                        y = df.venn.feces.metab_data, 
                                        use = "pairwise",method="spearman",adjust="BH", 
                                        alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(df.corr.feces.microb.metab.sp$r,   
         p.mat = df.corr.feces.microb.metab.sp$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')

r<- df.corr.feces.microb.metab.sp$r
p <- df.corr.feces.microb.metab.sp$p %>% as.data.frame()
cl <- df.sig.sp.DESeq2.MST %>%
  dplyr::rename(Species = feature) %>%
  dplyr::left_join(taxatable, by = "Species") %>%
  distinct(Species, .keep_all =TRUE)


ha1 = HeatmapAnnotation(`Enriched Group` = df.venn.feces.metab.select$GroupContrib, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(Class = cl$Class, show_annotation_name = F,
                    `Enriched Group` = df.sig.sp.DESeq2.MST$enrich_group, 
                    simple_anno_size = unit(0.3,"cm"),
                    col = list(`Enriched Group` = col1, Class = col2[unique(cl$Class)]))

p.feces.microb.metab.sp <- Heatmap(r, name = "Spearman's \u03C1",
                                show_heatmap_legend = TRUE,
                                top_annotation = ha1, left_annotation = ha2, 
                                cluster_rows = F, cluster_columns = F,
                                row_names_side = "left", column_names_side = "top",
                                row_names_gp = gpar(fontface = "italic"),
                                row_title = "Feces", row_title_side = "right",
                                row_title_rot = 90,
                                column_split = metab_level_group,
                                cluster_column_slices =F, column_title=NULL,
                                cell_fun = cell_fun, 
                                col = col_fun(9), column_names_rot = 45, 
                                width = ncol(r)*unit(6, "mm"), 
                                height = nrow(r)*unit(6, "mm")) 

####correlation between clinical and fecal metab####
df.corr.clinic.num.feces.metab <- corr.test(x = df.clinical.num.select.perio,
                                                y = df.venn.feces.metab_data, 
                                                use = "pairwise",method="spearman",adjust="BH", 
                                                alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corrplot(df.corr.clinic.num.feces.metab$r,   
         p.mat = df.corr.clinic.num.feces.metab$p.adj, tl.col = 'black', tl.srt = 45, 
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, na.label=" ",
         insig = 'label_sig')
write_xlsx(df.corr.clinic.num.feces.metab$p.adj %>% as.data.frame()%>% tibble::rownames_to_column(), 
           "Figures/Correlation/FD0.05/rel/df.corr.clinic.num.feces.metab.p.adj.xlsx")

r<- df.corr.clinic.num.feces.metab$r
p <- df.corr.clinic.num.feces.metab$p.adj %>% as.data.frame()



ha1 = HeatmapAnnotation(`Enriched Group` = df.venn.feces.metab.select$GroupContrib, 
                        show_annotation_name = F, show_legend = T,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))
ha2 = rowAnnotation(`Enriched Group` = df.clinical.num.group$Group, 
                        show_annotation_name = F, show_legend = F,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Enriched Group` = col1))

p.clinic.num.feces.metab <- Heatmap(r, name = "Spearman's \u03C1",
                                show_heatmap_legend = TRUE,
                                top_annotation = ha1, left_annotation = ha2, 
                                cluster_rows = F, cluster_columns = F,
                                row_names_side = "left", column_names_side = "top",
                                row_title = "", row_title_side = "right",
                                row_title_rot = 90,
                                column_split = metab_level_group,
                                cluster_column_slices =F, column_title=NULL,
                                cell_fun = cell_fun, 
                                col = col_fun(9), column_names_rot = 45, 
                                width = ncol(r)*unit(6, "mm"), 
                                height = nrow(r)*unit(6, "mm")) 

#### Boxplot of intersected genera of DESeq2 with and without BMI adjustment####
DESeq2.intersect.SBP.g <- intersect(sig.genus_DESeq2_SBP$feature,sig.genus_DESeq2_SBP.BMI$feature)
DESeq2.intersect.MSA.g <- intersect(sig.genus_DESeq2_MSA$feature,sig.genus_DESeq2_MSA.BMI$feature)
DESeq2.intersect.MST.g <- intersect(sig.genus_DESeq2_MST$feature,sig.genus_DESeq2_MST.BMI$feature)

df.log10.g.SBP.sel <- merge(sample_meta[sample_meta$Sample == "SBP",1:3], 
                            t(assay(tse.subsamples.g[[1]], "log10")[DESeq2.intersect.SBP.g,]),
                            by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names") 
df.log10.g.MSA.sel <- merge(sample_meta[sample_meta$Sample == "MSA",1:3], 
                            t(assay(tse.subsamples.g[[2]], "log10"))[,DESeq2.intersect.MSA.g],
                            by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names")
colnames(df.log10.g.MSA.sel)[4] <- DESeq2.intersect.MSA.g #only one MSA feature
df.log10.g.MST.sel <- merge(sample_meta[sample_meta$Sample == "MST",1:3], 
                                  t(assay(tse.subsamples.g[[3]], "log10")[DESeq2.intersect.MST.g,]),
                                  by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names")

plots.log10.g.SBP <- ggboxplot(melt(df.log10.g.SBP.sel), 
                               title = "Subgingival plaque",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, 
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,6)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 

plots.log10.g.MSA <- ggboxplot(melt(df.log10.g.MSA.sel), 
                               title = "Saliva",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, 
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,5)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 


plots.log10.g.MST <- ggboxplot(melt(df.log10.g.MST.sel), 
                               title = "Feces",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, 
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,6)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 

#horizonal plot
DESeq2.SBP.g <- sig.genus_DESeq2_SBP %>% arrange(ef_logFC)
DESeq2.MSA.g <- sig.genus_DESeq2_MSA %>% arrange(ef_logFC)
DESeq2.MST.g <- sig.genus_DESeq2_MST %>% arrange(ef_logFC)

df.log10.g.SBP.sel <- merge(sample_meta[sample_meta$Sample == "SBP",1:3], 
                            t(assay(tse.subsamples.g[[1]], "log10")[DESeq2.SBP.g$feature,]),
                            by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names") 
df.log10.g.MSA.sel <- merge(sample_meta[sample_meta$Sample == "MSA",1:3], 
                            t(assay(tse.subsamples.g[[2]], "log10")[DESeq2.MSA.g$feature,]),
                            by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names")

df.log10.g.MST.sel <- merge(sample_meta[sample_meta$Sample == "MST",1:3], 
                            t(assay(tse.subsamples.g[[3]], "log10")[DESeq2.MST.g$feature,]),
                            by = 'row.names', all = TRUE) %>%
  tibble::column_to_rownames("Row.names")

plots.log10.g.SBP <- ggboxplot(melt(df.log10.g.SBP.sel), 
                               title = "Subgingival plaque",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, orientation = "horizontal",
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,6)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 

plots.log10.g.MSA <- ggboxplot(melt(df.log10.g.MSA.sel), 
                               title = "Saliva",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, orientation = "horizontal",
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,5)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 


plots.log10.g.MST <- ggboxplot(melt(df.log10.g.MST.sel), 
                               title = "Feces",
                               x = "variable", y = "value",
                               color = "Disease", fill = "Disease", alpha = 0.2,
                               notch = FALSE, orientation = "horizontal",
                               add = "jitter") +
  labs(x ="", y="log10(Counts)") + 
  ylim(-1,6)+
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(size = 14, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 14),
        axis.title.y=element_text(size=14)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.signif", 
                     hide.ns = T, size =8, label.y = 4.5) 


# library(ggh4x)
# plots.log10.g.SBP.fix <- plots.log10.g.SBP + 
#                               force_panelsizes(rows = unit(8, "cm"),
#                                                                  cols = unit(24, "cm")) +
#                               theme(legend.position = "none")
# plots.log10.g.MSA.fix <- plots.log10.g.MSA + 
#                               force_panelsizes(rows = unit(8, "cm"),
#                                                cols = unit(1.5, "cm")) +
#                               theme(legend.position = "none")
# plots.log10.g.MST.fix <- plots.log10.g.MST + 
#                               force_panelsizes(rows = unit(8, "cm"),
#                                                cols = unit(12, "cm")) +
#                               theme(legend.position = "none")
# 
# library(gridExtra)
# lay <- rbind(c(1,1,1,1,1,1),
#              c(2,3,3))
# grid.arrange(plots.log10.g.SBP.fix,
#              plots.log10.g.MSA.fix, 
#              plots.log10.g.MST.fix, layout_matrix = lay) #need adjustment by Illustrator

egg::ggarrange(plots.log10.g.SBP + theme(legend.position = "none"), 
               plots.log10.g.MSA + theme(legend.position = "none"),
               plots.log10.g.MST + theme(legend.position = "none"),
               ncol = 1, heights = c(18, 2, 8))


#### microb correlation with picrust2 sig res####
p2.samples.list <- readRDS("PICRUSt2/picrust2_out_pipeline0.05/analysis/p2.samples.list.rds")
p2.samples.ANCOMBC2.res <- readRDS("PICRUSt2/picrust2_out_pipeline0.05/analysis/p2.samples.ANCOMBC2.res.rds")
df.SBP.g.select.clr <- t(as.data.frame(assay(tse.subsamples.g[[1]], "clr")))[,df.sig.genus.DESeq2.SBP$feature]
df.MSA.g.select.clr <- t(as.data.frame(assay(tse.subsamples.g[[2]], "clr")))[,df.sig.genus.DESeq2.MSA$feature]
df.MST.g.select.clr <- t(as.data.frame(assay(tse.subsamples.g[[3]], "clr")))[,df.sig.genus.DESeq2.MST$feature]


df.microb.g.select.counts.list <- list(SBP = df.SBP.g.select.clr,
                                       MSA = df.MSA.g.select.clr,
                                       MST = df.MST.g.select.clr) 
df.corr.microb.p2.summary.list <- list()
for (n in c("SBP", "MSA", "MST")){
  for (m in c("p2EC", "p2KO")){
    p2.res.sig <- p2.samples.ANCOMBC2.res[[n]][[paste0("tse.",m)]][["res_prim.sig"]] %>%
      mutate(ID = paste0(`function`, ": ", description)) 
    df.p2 <- p2.samples.list[[n]][[m]]
    df.p2.res.sig <- df.p2[p2.res.sig$`function`,] %>%
      tibble::rownames_to_column("function") %>%
      dplyr::left_join(p2.res.sig[,c("function", "ID")], by = "function") %>%
      tibble::column_to_rownames("ID") %>%
      dplyr::select(-`function`) %>% t() %>% as.data.frame()
    df.corr.microb.p2 <- corr.test(x = df.microb.g.select.counts.list[[n]],
                                   y = df.p2.res.sig,
                                   use = "pairwise",method="spearman",adjust="BH",
                                   alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    df.corr.microb.p2.r <- df.corr.microb.p2[["r"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
    colnames(df.corr.microb.p2.r) <- c("Genus", m, "Spearman_rho")
    df.corr.microb.p2.p <- df.corr.microb.p2[["p"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
    colnames(df.corr.microb.p2.p) <- c("Genus", m, "p_val")
    df.corr.microb.p2.p.adj <- df.corr.microb.p2[["p.adj"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
    colnames(df.corr.microb.p2.p.adj) <- c("Genus", m, "p_adj")
    df.corr.microb.p2.summary<-cbind(df.corr.microb.p2.r, df.corr.microb.p2.p, df.corr.microb.p2.p.adj)
    df.corr.microb.p2.summary.list[[n]][[m]] <- df.corr.microb.p2.summary[, c(1:3,6,9)]
    file.sig <- paste0("PICRUSt2/picrust2_out_pipeline0.05/analysis FD0.05/", 
                       paste(n,m, "corr.microb.p2.summary.csv", sep = "_"))
    write.csv(df.corr.microb.p2.summary.list[[n]][[m]], file.sig, row.names = F)
  }
}


for (n in c("SBP",  "MST")){
  p2.res.sig <- p2.samples.ANCOMBC2.res[[n]][["tse.p2PW"]][["res_prim.sig"]] %>%
    mutate(ID = paste0(pathway, ": ", description)) 
  df.p2 <- p2.samples.list[[n]][["p2PW"]]
  df.p2.res.sig <- df.p2[p2.res.sig$pathway,] %>%
    tibble::rownames_to_column("pathway") %>%
    dplyr::left_join(p2.res.sig[,c("pathway", "ID")], by = "pathway") %>%
    tibble::column_to_rownames("ID") %>%
    dplyr::select(-`pathway`) %>% t() %>% as.data.frame()
  df.corr.microb.p2 <- corr.test(x = df.microb.g.select.counts.list[[n]],
                                 y = df.p2.res.sig,
                                 use = "pairwise",method="spearman",adjust="BH",
                                 alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
  df.corr.microb.p2.r <- df.corr.microb.p2[["r"]] %>% as.data.frame %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(-rowname)
  colnames(df.corr.microb.p2.r) <- c("Genus", "p2PW", "Spearman_rho")
  df.corr.microb.p2.p <- df.corr.microb.p2[["p"]] %>% as.data.frame %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(-rowname)
  colnames(df.corr.microb.p2.p) <- c("Genus", "p2PW", "p_val")
  df.corr.microb.p2.p.adj <- df.corr.microb.p2[["p.adj"]] %>% as.data.frame %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(-rowname)
  colnames(df.corr.microb.p2.p.adj) <- c("Genus", "p2PW", "p_adj")
  df.corr.microb.p2.summary<-cbind(df.corr.microb.p2.r, df.corr.microb.p2.p, df.corr.microb.p2.p.adj)
  df.corr.microb.p2.summary.list[[n]][["p2PW"]] <- df.corr.microb.p2.summary[, c(1:3,6,9)]
  file.sig <- paste0("PICRUSt2/picrust2_out_pipeline0.05/analysis FD0.05/", 
                     paste(n,"p2PW", "corr.microb.p2.summary.csv", sep = "_"))
  write.csv(df.corr.microb.p2.summary.list[[n]][["p2PW"]], file.sig, row.names = F)
}


df.feces.g.ANCOMBC2 <- cbind(df.MST.g.select.clr,
                             df.p2.res.sig)
p_Co_P163_PWY <-ggscatter(df.feces.g.ANCOMBC2, x = "Coprococcus", 
                      y = "P163-PWY: L-lysine fermentation to acetate and butanoate", 
                      color = "black", shape = 21, size = 3, # Points color, shape and size
                      add = "reg.line", 
                      add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                      conf.int = TRUE, # Add confidence interval
                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                      cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                      xlab = "Coprococcus (clr)", ylab = "P163-PWY\n(functional aboundance)") 
p_Co_PWY_5677 <-ggscatter(df.feces.g.ANCOMBC2, x = "Coprococcus", 
                          y = "PWY-5677: succinate fermentation to butanoate", 
                          color = "black", shape = 21, size = 3, # Points color, shape and size
                          add = "reg.line", 
                          add.params = list(color = "#00468B", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE, # Add confidence interval
                          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                          xlab = "Coprococcus (clr)", ylab = "PWY-5677\n(functional aboundance)") 

p_Akk_P163_PWY<- ggscatter(df.feces.g.ANCOMBC2, x = "Akkermansia", 
                           y = "P163-PWY: L-lysine fermentation to acetate and butanoate", 
                           color = "black", shape = 21, size = 3, # Points color, shape and size
                           add = "reg.line", 
                           add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                           conf.int = TRUE, # Add confidence interval
                           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                           cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                           xlab = "Akkermansia (clr)", ylab = "P163-PWY\n(functional aboundance)") 
  
p_Akk_PWY_5677 <- ggscatter(df.feces.g.ANCOMBC2, x = "Akkermansia", 
                          y = "PWY-5677: succinate fermentation to butanoate", 
                          color = "black", shape = 21, size = 3, # Points color, shape and size
                          add = "reg.line", 
                          add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE, # Add confidence interval
                          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                          xlab = "Akkermansia (clr)", ylab = "PWY-5677\n(functional aboundance)") 

p_La_P163_PWY <- ggscatter(df.feces.g.ANCOMBC2, x = "Lachnoclostridium", 
                           y = "P163-PWY: L-lysine fermentation to acetate and butanoate", 
                           color = "black", shape = 21, size = 3, # Points color, shape and size
                           add = "reg.line", 
                           add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                           conf.int = TRUE, # Add confidence interval
                           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                           cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                           xlab = "Lachnoclostridium (clr)", ylab = "P163-PWY\n(functional aboundance)") 
  
p_La_PWY_5677 <- ggscatter(df.feces.g.ANCOMBC2, x = "Lachnoclostridium", 
                           y = "PWY-5677: succinate fermentation to butanoate", 
                           color = "black", shape = 21, size = 3, # Points color, shape and size
                           add = "reg.line", 
                           add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                           conf.int = TRUE, # Add confidence interval
                           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                           cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
                           xlab = "Lachnoclostridium (clr)", ylab = "PWY-5677\n(functional aboundance)")
ggarrange(p_Co_P163_PWY, p_Co_PWY_5677,
          p_Akk_P163_PWY, p_Akk_PWY_5677, ncol = 4, nrow = 1)
