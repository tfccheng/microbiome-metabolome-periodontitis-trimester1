library(ANCOMBC)
library(mia)
library(dplyr)
library(tidyr)
#library(knitr)
library(DT)

set.seed(123)
setwd("")
tse.micro <- readRDS("") #read processed data tse file


tse.micro.subsamples <- list()
for (n in c("SBP", "MSA", "MST")){
  tse.micro.subsample <- tse.micro[, tse.micro$Sample %in% n]
  tse.micro.subsamples[[n]] <- tse.micro.subsample
}

output.subsamples <- list()
for (n in c("SBP", "MSA", "MST")){
  for (m in c("Family", "Genus", "Species")){
    output = ancombc2(data = tse.micro.subsamples[[n]], assay_name = "counts", tax_level = m,
                      fix_formula = "Disease", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      group = "Disease", struc_zero = TRUE, neg_lb = TRUE,
                      alpha = 0.2, n_cl = 12, verbose = TRUE)
    output.subsamples[[n]][[m]] <- output
  }
}
for (n in c("SBP", "MSA", "MST")){
  for (m in c("Family", "Genus", "Species")){
    write.csv(output.subsamples[[n]][[m]]$res, 
              paste0("Figures/ANCOM-BC2/ANCOM-BC2_output_", n, "_", m, ".csv"), row.names = F)
    
  }
}



output.subsamples.BMI <- list()
for (n in c("SBP", "MSA", "MST")){
  for (m in c("Family", "Genus", "Species")){
    output = ancombc2(data = tse.micro.subsamples[[n]], assay_name = "counts", tax_level = m,
                      fix_formula = "Disease + BMI_T1.x", rand_formula = NULL,
                      p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      group = "Disease", struc_zero = TRUE, neg_lb = TRUE,
                      alpha = 0.2, n_cl = 12, verbose = TRUE)
    output.subsamples.BMI[[n]][[m]] <- output
  }
}

for (n in c("SBP", "MSA", "MST")){
  for (m in c("Family", "Genus", "Species")){
    write.csv(output.subsamples.BMI[[n]][[m]]$res, 
              paste0("Figures/ANCOM-BC2/ANCOM-BC2_output_BMI_", n, "_", m, ".csv"), row.names = F)
    
  }
}







  