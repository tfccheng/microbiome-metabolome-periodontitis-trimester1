library("data.table")   
library("phyloseq")
library("ALDEx2")
library("tidyverse")
library(readxl)
library(dplyr)
library(TreeSummarizedExperiment)
library(mia)
library(ANCOMBC)
library(KEGGREST)

setwd("./PICRUSt2")
mapfile = "description_mapfiles"
list.files(mapfile, recursive = TRUE)
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_KO_BRITE = paste0(mapfile, "/picrust1_KO_BRITE_map.tsv")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

p2_EC = "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
p2_KO = "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
p2_PW = "/pathways_out/path_abun_unstrat.tsv.gz"

mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapKO_BRITE = read.table(mapfile_KO_BRITE, header=TRUE, sep="\t", quote = "", 
                         stringsAsFactors = FALSE, comment.char="", row.names=1) %>%
  tibble::rownames_to_column("function")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")


p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC) %>% as.data.frame()

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO) %>% as.data.frame()

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW) %>% as.data.frame()

meta <- read.csv("meta new.csv",
                 row.names = 1)
meta <- meta[ order(row.names(meta)), ]
##by samples
#Subgingival (SBP), saliva (MSA) and feces (MST)
p2.samples.list <- list()
for (n in c("SBP", "MSA", "MST")){
  res <- subset(meta, Sample == n)
  df.p2EC <- p2EC[, which(names(p2EC) %in% row.names(res))]
  df.p2KO <- p2KO[, which(names(p2KO) %in% row.names(res))]
  df.p2PW <- p2PW[, which(names(p2PW) %in% row.names(res))]
  p2.samples.list[[n]]$p2EC <- df.p2EC
  p2.samples.list[[n]]$p2KO <- df.p2KO
  p2.samples.list[[n]]$p2PW <- df.p2PW
  p2.samples.list[[n]]$meta <- res
  assays.EC = S4Vectors::SimpleList(counts = as.matrix(df.p2EC))
  assays.KO = S4Vectors::SimpleList(counts = as.matrix(df.p2KO))
  assays.PW = S4Vectors::SimpleList(counts = as.matrix(df.p2PW))
  p2.samples.list[[n]][["tse.p2EC"]] <- TreeSummarizedExperiment(assays = assays.EC, colData = res) #tse object for ANCOMBC2
  p2.samples.list[[n]][["tse.p2KO"]] <- TreeSummarizedExperiment(assays = assays.KO, colData = res)
  p2.samples.list[[n]][["tse.p2PW"]] <- TreeSummarizedExperiment(assays = assays.PW, colData = res)
}
saveRDS(p2.samples.list, "analysis/p2.samples.list.rds")
set.seed(1234)

p2.samples.ANCOMBC2.list <- list()
for (n in c("SBP", "MSA", "MST")){
  for (m in c("tse.p2EC", "tse.p2KO", "tse.p2PW")){
    tse = p2.samples.list[[n]][m]
    output = ancombc2(data = p2.samples.list[[n]][[m]], assay_name = "counts",
                      fix_formula = "Disease", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      group = "Disease", struc_zero = TRUE, neg_lb = TRUE,
                      alpha = 0.2, n_cl = 12, verbose = TRUE)
    p2.samples.ANCOMBC2.list[[n]][[m]] <- output
  }
}


p2.samples.ANCOMBC2.res <- list()
for (n in c("SBP", "MSA", "MST")){
  for (m in c("tse.p2EC", "tse.p2KO", "tse.p2PW")){
    if (m == "tse.p2EC") {
      res_prim <- p2.samples.ANCOMBC2.list[[n]][[m]]$res %>%
        rename("taxon" = "function")
      res_prim <- dplyr::inner_join(res_prim, mapEC, by = "function") %>%
        arrange(p_DiseasePD)
    } else if (m == "tse.p2KO") {
      res_prim <- p2.samples.ANCOMBC2.list[[n]][[m]]$res %>%
        rename("taxon" = "function")
      res_prim <- dplyr::inner_join(res_prim, mapKO, by = "function") %>%
        arrange(p_DiseasePD)
    } else{
      res_prim <- p2.samples.ANCOMBC2.list[[n]][[m]]$res %>%
        rename("taxon" = "pathway")
      res_prim <- dplyr::inner_join(res_prim, mapPW, by = "pathway") %>%
        arrange(p_DiseasePD)
    }
    p2.samples.ANCOMBC2.res[[n]][[m]]$res_prim <- res_prim
    res_prim.sig <- subset(res_prim, p_DiseasePD < 0.05)
    p2.samples.ANCOMBC2.res[[n]][[m]]$res_prim.sig <- res_prim.sig
    file <- paste0("analysis/", 
                  paste(n,m, "ANCOMBC2_res.csv", sep = "_"))
    write.csv(res_prim[-(12:13)], file, col.names = T, row.names = F)
    file.sig <- paste0("analysis/", 
                      paste(n,m, "ANCOMBC2_res.sig.csv", sep = "_"))
    write.csv(res_prim.sig[-(12:13)], file.sig, col.names = T, row.names = F)
  }
}

saveRDS(p2.samples.ANCOMBC2.res, "analysis/p2.samples.ANCOMBC2.res.rds")


for (n in c("SBP", "MSA", "MST")){
    res_prim.sig <- p2.samples.ANCOMBC2.res[[n]][["tse.p2KO"]]$res_prim.sig
    a<- dplyr::left_join(res_prim.sig, mapKO_BRITE[,c(1,3)], by = "function")
    file.sig <- paste0("analysis/", 
                      paste(n,"tse.p2KO", "ANCOMBC2_res.sig.full_KEGG.csv", sep = "_"))
    write.csv(a[-(12:13)], file.sig, col.names = T, row.names = F)
}







