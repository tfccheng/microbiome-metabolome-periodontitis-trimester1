library(dplyr)
library(readxl)
library(mia)
library(stringr)
library(ggplot2)
library(igraph)

setwd("")

####metabolome data preparation####
# load IQR filtered metabolome
df_metab_feces_X.log.pareto <- read.csv("Data filtering/feces/data_normalized_no_label.csv") %>%
  tibble::column_to_rownames("X") %>% t() %>% as.data.frame()
df_metab_saliva_X.log.pareto <- read.csv("Data filtering/saliva/data_normalized_no_label.csv") %>%
  tibble::column_to_rownames("X") %>% t() %>% as.data.frame()
df_metab_serum_X.log.pareto <- read.csv("Data filtering/serum/data_normalized_no_label.csv") %>%
  tibble::column_to_rownames("X") %>% t() %>% as.data.frame()

df.ttest.feces <- read.csv("Data filtering/feces/t_test.csv") %>%
  tibble::column_to_rownames("X")
df.ttest.saliva <- read.csv("Data filtering/saliva/t_test.csv") %>%
  tibble::column_to_rownames("X")
df.ttest.serum <- read.csv("Data filtering/serum/t_test.csv") %>%
  tibble::column_to_rownames("X")
df.wilcox.feces <- read.csv("Data filtering/feces/wilcox_rank.csv") %>%
  tibble::column_to_rownames("X")
df.wilcox.saliva <- read.csv("Data filtering/saliva/wilcox_rank.csv") %>%
  tibble::column_to_rownames("X")
df.wilcox.serum <- read.csv("Data filtering/serum/wilcox_rank.csv") %>%
  tibble::column_to_rownames("X")

df.ttest.feces.all <- read.csv("Data filtering/feces/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)
df.ttest.saliva.all <- read.csv("Data filtering/saliva/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)
df.ttest.serum.all <- read.csv("Data filtering/serum/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)


df_metab_feces_Case_Ctrl.all <- read_xlsx("Data filtering/feces/Case_Ctrl.all.xlsx") %>%
  filter(ID %in% df.ttest.feces.all$ID) 
df_metab_saliva_Case_Ctrl.all <- read_xlsx("Data filtering/saliva/Case_Ctrl.all.xlsx") %>%
  filter(ID %in% df.ttest.saliva.all$ID) 
df_metab_serum_Case_Ctrl.all <- read_xlsx("Data filtering/serum/Case_Ctrl.all.xlsx") %>%
  filter(ID %in% df.ttest.serum.all$ID) 

#prepare for functional analysis by MetaboAnalyst
df_metab_feces <- df.ttest.feces.all %>%
  dplyr::left_join(df_metab_feces_Case_Ctrl.all, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT)
write.csv(df_metab_feces, "Data filtering/feces/df_metab_feces_for_func.csv", row.names = F)
df_metab_saliva <- df.ttest.saliva.all %>%
  dplyr::left_join(df_metab_saliva_Case_Ctrl.all, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT)
write.csv(df_metab_saliva, "Data filtering/saliva/df_metab_saliva_for_func.csv", row.names = F)
df_metab_serum <- df.ttest.serum.all %>%
  dplyr::left_join(df_metab_serum_Case_Ctrl.all, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT)
write.csv(df_metab_serum, "Data filtering/serum/df_metab_serum_for_func.csv", row.names = F)


#filter most sig between ttest and wilcox
df_metab_feces_X.log.pareto.sig <- df_metab_feces_X.log.pareto[,rownames(df.ttest.feces)]
df_metab_saliva_X.log.pareto.sig <- df_metab_saliva_X.log.pareto[,rownames(df.ttest.saliva)]
df_metab_serum_X.log.pareto.sig <- df_metab_serum_X.log.pareto[,rownames(df.wilcox.serum)]


#identification
df_metab_feces_identification <- read_xlsx("Metab identification/feces/identification.xlsx") %>%
  na_if("-")
df_metab_feces_identification$ident <- ifelse(is.na(df_metab_feces_identification$MS2Metabolite), 
                                              ifelse(is.na(df_metab_feces_identification$MS1hmdbName), 
                                                     ifelse(is.na(df_metab_feces_identification$MS1keggName), NA, 
                                                            df_metab_feces_identification$MS1keggName),
                                                     df_metab_feces_identification$MS1hmdbName),
                                              df_metab_feces_identification$MS2Metabolite)
df_metab_feces_identification$identification <- ifelse(is.na(df_metab_feces_identification$MS2Metabolite), 
                                                       ifelse(is.na(df_metab_feces_identification$MS1hmdbName), 
                                                              ifelse(is.na(df_metab_feces_identification$MS1keggName), NA, 
                                                                     paste(df_metab_feces_identification$MS1keggName, " (MS1 ", 
                                                                           df_metab_feces_identification$MS1keggID,")", sep = "")),
                                                              paste(df_metab_feces_identification$MS1hmdbName, " (MS1 ",
                                                                    df_metab_feces_identification$MS1hmdbID,")", sep = "")),
                                                       paste(df_metab_feces_identification$MS2Metabolite, " (MS2 ",
                                                             df_metab_feces_identification$MS2hmdb,")", sep = ""))
df_metab_saliva_identification <- read_xlsx("Metab identification/saliva/identification.xlsx") %>%
  na_if("-")
df_metab_saliva_identification$ident <- ifelse(is.na(df_metab_saliva_identification$MS2Metabolite), 
                                               ifelse(is.na(df_metab_saliva_identification$MS1hmdbName), 
                                                      ifelse(is.na(df_metab_saliva_identification$MS1keggName), NA, 
                                                             df_metab_saliva_identification$MS1keggName),
                                                      df_metab_saliva_identification$MS1hmdbName),
                                               df_metab_saliva_identification$MS2Metabolite)
df_metab_saliva_identification$identification <- ifelse(is.na(df_metab_saliva_identification$MS2Metabolite), 
                                                        ifelse(is.na(df_metab_saliva_identification$MS1hmdbName), 
                                                               ifelse(is.na(df_metab_saliva_identification$MS1keggName), NA, 
                                                                      paste(df_metab_saliva_identification$MS1keggName, " (MS1 ", 
                                                                            df_metab_saliva_identification$MS1keggID,")", sep = "")),
                                                               paste(df_metab_saliva_identification$MS1hmdbName, " (MS1 ",
                                                                     df_metab_saliva_identification$MS1hmdbID,")", sep = "")),
                                                        paste(df_metab_saliva_identification$MS2Metabolite, " (MS2 ",
                                                              df_metab_saliva_identification$MS2hmdb,")", sep = ""))
df_metab_serum_identification <- read_xlsx("Metab identification/serum/identification.xlsx") %>%
  na_if("-")
df_metab_serum_identification$ident <- ifelse(is.na(df_metab_serum_identification$MS2Metabolite), 
                                              ifelse(is.na(df_metab_serum_identification$MS1hmdbName), 
                                                     ifelse(is.na(df_metab_serum_identification$MS1keggName), NA, 
                                                            df_metab_serum_identification$MS1keggName),
                                                     df_metab_serum_identification$MS1hmdbName),
                                              df_metab_serum_identification$MS2Metabolite)
df_metab_serum_identification$identification <- ifelse(is.na(df_metab_serum_identification$MS2Metabolite), 
                                                       ifelse(is.na(df_metab_serum_identification$MS1hmdbName), 
                                                              ifelse(is.na(df_metab_serum_identification$MS1keggName), NA, 
                                                                     paste(df_metab_serum_identification$MS1keggName, " (MS1 ", 
                                                                           df_metab_serum_identification$MS1keggID,")", sep = "")),
                                                              paste(df_metab_serum_identification$MS1hmdbName, " (MS1 ",
                                                                    df_metab_serum_identification$MS1hmdbID,")", sep = "")),
                                                       paste(df_metab_serum_identification$MS2Metabolite, " (MS2 ",
                                                             df_metab_serum_identification$MS2hmdb,")", sep = ""))
saveRDS(df_metab_feces_identification, "Data prepare/df_metab_feces_identification.rds")
saveRDS(df_metab_saliva_identification, "Data prepare/df_metab_saliva_identification.rds")
saveRDS(df_metab_serum_identification, "Data prepare/df_metab_serum_identification.rds")



dim(df_metab_feces_X.log.pareto)
metab.meta <- read_excel("metabolome sample meta.xlsx")

####microbiome data preparation####
## microbiome data of sugbgingival plaque, saliva and feces
tse.micro <- readRDS("D:/Research Data/SZ project/T1-Result-X101SC20072895-Z01-F001-B1-10/personal analysis/T1 microbiome analysis/tse.silva.wclean.clinical.calprotectin.rds")

write.csv(data.frame(rowData(tse.micro)), "Data prepare/taxa.micro.csv")
write.csv(data.frame(rowData(agglomerateByRank(tse.micro, rank = "Genus", onRankOnly = T))), "Data prepare/taxa.micro.genus.csv")
write.csv(data.frame(rowData(agglomerateByRank(tse.micro, rank = "Species", onRankOnly = T))), "Data prepare/taxa.micro.species.csv")
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

df.micro.saliva <- as.data.frame(assay(tse.micro.subsamples[["MSA"]], "rclr")) %>% t() %>% as.data.frame()
df.micro.feces <- as.data.frame(assay(tse.micro.subsamples[["MST"]], "rclr")) %>% t() %>% as.data.frame()
df.micro.plaque <- as.data.frame(assay(tse.micro.subsamples[["SBP"]], "rclr")) %>% t() %>% as.data.frame()


df.micro.saliva.g <- altExp(tse.micro.subsamples[["MSA"]], "Genus") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()
df.micro.feces.g <- altExp(tse.micro.subsamples[["MST"]], "Genus") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()
df.micro.plaque.g <- altExp(tse.micro.subsamples[["SBP"]], "Genus") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()

df.micro.saliva.s <- altExp(tse.micro.subsamples[["MSA"]], "Species") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()
df.micro.feces.s <- altExp(tse.micro.subsamples[["MST"]], "Species") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()
df.micro.plaque.s <- altExp(tse.micro.subsamples[["SBP"]], "Species") %>%
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% 
  assay("rclr") %>%
  as.data.frame() %>% t() %>% as.data.frame()

micro.meta.saliva <- as.data.frame(colData(tse.micro.subsamples[["MSA"]]))
micro.meta.feces <- as.data.frame(colData(tse.micro.subsamples[["MST"]]))
micro.meta.plaque <- as.data.frame(colData(tse.micro.subsamples[["SBP"]]))

####clinical data preparation####
## clinical continuous data
df.clinical <- as.data.frame(colData(tse.micro.subsamples[["MSA"]]))[,-c(1:3,5)] %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("ID")
df.clinical.num <- df.clinical[,54:80]
#impute NA
nipals.tune.clinical = nipals(df.clinical.num, ncomp = 10)$eig
barplot(nipals.tune.clinical, xlab = 'Principal component', ylab = 'Explained variance')
df.clinical.num.impute <- impute.nipals(df.clinical.num, ncomp = 10) # impute NAs if only several NAs for each item
#remove BMI, SBP, DBP
df.clinical.num.impute.2 <- df.clinical.num.impute[,-(1:3)]

meta.metab.micro.combine.feces <- dplyr::inner_join(micro.meta.feces[1:5], metab.meta, by = "ID")
meta.metab.micro.combine.saliva <- dplyr::inner_join(micro.meta.saliva[1:5], metab.meta, by = "ID")
saveRDS(meta.metab.micro.combine.feces, "Data prepare/meta.metab.micro.combine.feces.rds")
saveRDS(meta.metab.micro.combine.saliva, "Data prepare/meta.metab.micro.combine.saliva.rds")
#for metaboanalyst multivariate analysis meta
meta_clinic_for_metaboanalyst <- df.clinical.num.impute %>% tibble::rownames_to_column("clinicaldataID") %>%
  dplyr::left_join(meta.metab.micro.combine.feces, by = "clinicaldataID") %>%
  dplyr::select(-Sample, -clinicaldataID, -order, -class, -ID, -Source) %>% 
  dplyr::rename(Sample = sample_ID) %>%
  relocate(Sample, .before = 1) %>% relocate(Disease, .after = 1)
write.csv(meta_clinic_for_metaboanalyst, "meta_clinic_for_metaboanalyst.csv", row.names = F)

####sort samples####
df_metab_feces_X.log.pareto <- df_metab_feces_X.log.pareto[meta.metab.micro.combine.feces$sample_ID,]
identical(rownames(df_metab_feces_X.log.pareto), meta.metab.micro.combine.feces$sample_ID) # check sample match
rownames(df_metab_feces_X.log.pareto) <- meta.metab.micro.combine.feces$clinicaldataID 

df_metab_saliva_X.log.pareto <- df_metab_saliva_X.log.pareto[meta.metab.micro.combine.saliva$sample_ID,]
identical(rownames(df_metab_saliva_X.log.pareto), meta.metab.micro.combine.saliva$sample_ID)
rownames(df_metab_saliva_X.log.pareto) <- meta.metab.micro.combine.saliva$clinicaldataID 

df_metab_serum_X.log.pareto <- df_metab_serum_X.log.pareto[meta.metab.micro.combine.saliva$sample_ID,]
identical(rownames(df_metab_serum_X.log.pareto), meta.metab.micro.combine.saliva$sample_ID)
rownames(df_metab_serum_X.log.pareto) <- meta.metab.micro.combine.saliva$clinicaldataID 
#sig
df_metab_feces_X.log.pareto.sig <- df_metab_feces_X.log.pareto.sig[meta.metab.micro.combine.feces$sample_ID,]
identical(rownames(df_metab_feces_X.log.pareto.sig), meta.metab.micro.combine.feces$sample_ID) # check sample match
rownames(df_metab_feces_X.log.pareto.sig) <- meta.metab.micro.combine.feces$clinicaldataID 

df_metab_saliva_X.log.pareto.sig <- df_metab_saliva_X.log.pareto.sig[meta.metab.micro.combine.saliva$sample_ID,]
identical(rownames(df_metab_saliva_X.log.pareto.sig), meta.metab.micro.combine.saliva$sample_ID)
rownames(df_metab_saliva_X.log.pareto.sig) <- meta.metab.micro.combine.saliva$clinicaldataID 

df_metab_serum_X.log.pareto.sig <- df_metab_serum_X.log.pareto.sig[meta.metab.micro.combine.saliva$sample_ID,]
identical(rownames(df_metab_serum_X.log.pareto.sig), meta.metab.micro.combine.saliva$sample_ID)
rownames(df_metab_serum_X.log.pareto.sig) <- meta.metab.micro.combine.saliva$clinicaldataID 

saveRDS(df_metab_feces_X.log.pareto, "Data prepare/df_metab_feces_X.log.pareto.rds")
saveRDS(df_metab_saliva_X.log.pareto, "Data prepare/df_metab_saliva_X.log.pareto.rds")
saveRDS(df_metab_serum_X.log.pareto, "Data prepare/df_metab_serum_X.log.pareto.rds")

saveRDS(df_metab_feces_X.log.pareto.sig, "Data prepare/df_metab_feces_X.log.pareto.sig.rds")
saveRDS(df_metab_saliva_X.log.pareto.sig, "Data prepare/df_metab_saliva_X.log.pareto.sig.rds")
saveRDS(df_metab_serum_X.log.pareto.sig, "Data prepare/df_metab_serum_X.log.pareto.sig.rds")
df.clinical.num.impute <- df.clinical.num.impute[meta.metab.micro.combine.saliva$ID,]
rownames(df.clinical.num.impute) <- meta.metab.micro.combine.saliva$clinicaldataID 
colnames(df.clinical.num.impute) <- c("BMI", "SBP", "DBP", "WBC", "CRP",        
                                      "MA", "UCr", "MA/UCr", "Uric acid", "TBIL", 
                                      "TBA", "ALT", "AST", "TP", "ALB",         
                                      "TC", "TG", "HDLc", "LDLc", "Ferritin",         
                                      "Insulin", "Blood glucose", "HbA1c", "TSH", "FT4", "TPO-Ab","Calprotectin")
df.clinical.num.impute.2 <- df.clinical.num.impute.2[meta.metab.micro.combine.saliva$ID,]
rownames(df.clinical.num.impute.2) <- meta.metab.micro.combine.saliva$clinicaldataID 
colnames(df.clinical.num.impute.2) <- c("WBC", "CRP",        
                                        "MA", "UCr", "MA/UCr", "Uric acid", "TBIL", 
                                        "TBA", "ALT", "AST", "TP", "ALB",         
                                        "TC", "TG", "HDLc", "LDLc", "Ferritin",         
                                        "Insulin", "Blood glucose", "HbA1c", "TSH", "FT4", "TPO-Ab", "Calprotectin")
saveRDS(df.clinical.num.impute, "Data prepare/df.clinical.num.impute.rds")
saveRDS(df.clinical.num.impute.2, "Data prepare/df.clinical.num.impute.2.rds")

rownames(df.micro.saliva) <- meta.metab.micro.combine.saliva$clinicaldataID
rownames(df.micro.feces) <- meta.metab.micro.combine.feces$clinicaldataID
rownames(df.micro.plaque) <- meta.metab.micro.combine.feces$clinicaldataID
rownames(df.micro.saliva.g) <- meta.metab.micro.combine.saliva$clinicaldataID
rownames(df.micro.feces.g) <- meta.metab.micro.combine.feces$clinicaldataID
rownames(df.micro.plaque.g) <- meta.metab.micro.combine.feces$clinicaldataID
rownames(df.micro.saliva.s) <- meta.metab.micro.combine.saliva$clinicaldataID
rownames(df.micro.feces.s) <- meta.metab.micro.combine.feces$clinicaldataID
rownames(df.micro.plaque.s) <- meta.metab.micro.combine.feces$clinicaldataID

saveRDS(df.micro.saliva, "Data prepare/df.micro.saliva.rds")
saveRDS(df.micro.feces, "Data prepare/df.micro.feces.rds")
saveRDS(df.micro.plaque, "Data prepare/df.micro.plaque.rds")
saveRDS(df.micro.saliva.g, "Data prepare/df.micro.saliva.g.rds")
saveRDS(df.micro.feces.g, "Data prepare/df.micro.feces.g.rds")
saveRDS(df.micro.plaque.g, "Data prepare/df.micro.plaque.g.rds")
saveRDS(df.micro.saliva.s, "Data prepare/df.micro.saliva.s.rds")
saveRDS(df.micro.feces.s, "Data prepare/df.micro.feces.s.rds")
saveRDS(df.micro.plaque.s, "Data prepare/df.micro.plaque.s.rds")