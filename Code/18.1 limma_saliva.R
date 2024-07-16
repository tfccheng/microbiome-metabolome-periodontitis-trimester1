library(MetaboAnalystR)
setwd("")
mSet<-InitDataObjects("pktable", "mf", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, "data_normalized_no_label for metaboanalyst.saliva.csv", "colmf");
mSet<-ReadMetaData(mSet, "meta_clinic_for_metaboanalyst.csv");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckMeta(mSet, 1)
mSet<-SetDataTypeOfMeta(mSet)
FilterVariable(mSet, "none", 0, "F", F)
mSet<-SanityCheckData(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
adj.vec <<- "BMI"
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_0_dpi72.pdf", "pdf", "Disease", "NPD", "NA" , 0.05,"PD")
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_0_dpi72.pdf", "pdf", "Disease", "NPD", "NA" , 1,"PD") #all values

remove(adj.vec)
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_0_dpi72.pdf", "pdf", "Disease", "NPD", "NA" , 1,"PD") #all values without BMI adj




df.covariate.BMI <- read.csv("covariate_result.saliva.BMI.all.csv") %>%
  rename(p.BMI = P.Value)
df.covariate <- read.csv("covariate_result.saliva.all.csv") %>%
  rename(p = P.Value)

df.p.value <- df.covariate %>%
  merge(df.covariate.BMI, by = "X") %>%
  dplyr::select(X, p, p.BMI) %>%
  mutate(significance = ifelse(p < 0.05 & p.BMI < 0.05, "Significant", 
                               ifelse(p < 0.05 & p.BMI >= 0.05, "Non-sig. BMI adjusted", 
                                      ifelse(p >= 0.05 & p.BMI < 0.05, "Sig. BMI adjusted", "Non-sig.")))) %>%
  mutate(pval.no = -log10(p)) %>% mutate(pval.adj = -log10(p.BMI)) %>%
  arrange(p)

df.p.value$significance <- factor(df.p.value$significance, 
                                  levels = c("Significant", "Sig. BMI adjusted",
                                             "Non-sig.", "Non-sig. BMI adjusted"))

summary(factor(df.p.value$significance))


p<- ggplot(df.p.value, mapping = aes(x = pval.no, y = pval.adj, col = significance, label = X)) +
  guides(size="none") +
  geom_point(aes(size=pval.adj), alpha=0.5) +
  scale_color_manual(values=c("#0070ff", "#158315", "grey", "#c4002f"), name = "",
                     guide = guide_legend(override.aes = list(size = 3))) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey20") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20") +
  xlab("-log10(P-value): no covariate adjustment") +
  ylab("-log10(P-value): adjusted") +
  geom_text_repel(data = df.p.value[c(1:5),], 
                  aes(x=pval.no,y=pval.adj,label=X), show.legend = F) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "top")

df_metab_saliva_identification <- readRDS("Data filtering/df_metab_saliva_identification.rds")
df.covariate.BMI.sig <- read.csv("covariate_result.saliva.BMI.csv") %>%
  rename(p.BMI = P.Value)
df.covariate.BMI_w.anno <- df.covariate.BMI.sig %>% rename(ID = X) %>%
  dplyr::left_join(df_metab_saliva_identification, by = "ID") %>% 
  arrange(desc(IsMS2ident), desc(p.BMI))
write.csv(df.covariate.BMI_w.anno, "df.covariate.BMI_w.anno.saliva.csv", row.names = F)
pos<-df.covariate.BMI_w.anno %>% filter(IsMS2ident == 1 & logFC>0)
neg<-df.covariate.BMI_w.anno %>% filter(IsMS2ident == 1 & logFC<0)
summary(factor(pos$MS2superclass))
summary(factor(neg$MS2superclass))

df_metab_saliva_Case_Ctrl.all <- read_xlsx("Data filtering/saliva/Case_Ctrl.all.xlsx") %>%
  filter(ID %in% df.p.value$X) 
saveRDS(df_metab_saliva_Case_Ctrl.all%>% filter(MS2superclass != "Unknown") %>% dplyr::select(MS2superclass), 
        "saliva_MS2superclass.rds")

df_metab_saliva_Case_Ctrl.all.limma.sig <- df_metab_saliva_Case_Ctrl.all %>%
  filter(ID %in% df.covariate.BMI.sig$X) 

df.ttest.saliva.all <- read.csv("Data filtering/saliva/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)

#prepare for functional analysis by MetaboAnalyst
df_metab_saliva.limma.sig <- df.ttest.saliva.all %>%
  dplyr::right_join(df_metab_saliva_Case_Ctrl.all.limma.sig, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT) %>%
  arrange(mode)
write.csv(df_metab_saliva.limma.sig, "df_metab_saliva.limma.sig_for_func.csv", row.names = F)

#functional enrichment mummichog and GSEA
mSet_fun <-InitDataObjects("mass_all", "mummichog", FALSE)
mSet_fun<-SetPeakFormat(mSet_fun, "rmp")
mSet_fun<-UpdateInstrumentParameters(mSet_fun, 10.0, "mixed", "yes", 0.02);
mSet_fun<-Read.PeakListData(mSet_fun, "df_metab_saliva.limma.sig_for_func.csv");
mSet_fun<-SanityCheckMummichogData(mSet_fun)


mSet_fun<-SetPeakEnrichMethod(mSet_fun, "integ", "v2")
mSet_fun<-SetMummichogPval(mSet_fun, 0.01)
mSet_fun<-PerformPSEA(mSet_fun, "hsa_mfn", "current", 3 , 100)
mSet_fun<-PlotPSEAIntegPaths(mSet_fun, "func_enrich_saliva_mummichog2_GSEA", "pdf", 72, width=NA)
