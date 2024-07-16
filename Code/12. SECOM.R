library(MicrobiomeAnalystR)

mbSet<-Init.mbSetObj()
mbSet<-SetModuleType(mbSet, "mdp")
mbSet<-ReadSampleTable(mbSet, "meta new.csv");
mbSet<-Read16STaxaTable(mbSet, "ASV_silva_species_combined_genus.txt");
mbSet<-Read16SAbundData(mbSet, "ASVs_counts.txt","text","SILVA","T","false");
mbSet<-SanityCheckData(mbSet, "text");
mbSet<-SanityCheckSampleData(mbSet);
mbSet<-SetMetaAttributes(mbSet, "1")
mbSet<-PlotLibSizeView(mbSet, "norm_libsizes_0","png");
mbSet<-CreatePhyloseqObj(mbSet, "text","SILVA","F" , "false")
mbSet<-ApplyAbundanceFilter(mbSet, "prevalence", 4, 0.2);
mbSet<-ApplyVarianceFilter(mbSet, "iqr", 0.1);
smpl.nm.vec <- c("GCF.T1.01","GCF.T1.02","GCF.T1.03","GCF.T1.04","GCF.T1.05","GCF.T1.06","GCF.T1.09","GCF.T1.10","GCF.T1.11","GCF.T1.13","GCF.T1.14","GCF.T1.16","GCF.T1.17","GCF.T1.18","GCF.T1.19","GCF.T1.21","GCF.T1.22","GCF.T1.23","GCF.T1.24","GCF.T1.25","GCF.T1.26","GCF.T1.27","GCF.T1.28","GCF.T1.29","GCF.T1.30","GCF.T1.31","GCF.T1.32","GCF.T1.33","GCF.T1.34","GCF.T1.35","GCF.T1.36","GCF.T1.39","GCF.T1.40","GCF.T1.41","GCF.T1.42","GCF.T1.44","GCF.T1.45","GCF.T1.46","GCF.T1.47","GCF.T1.49","GCF.T1.50","GCF.T1.51","GCF.T1.52","GCF.T1.55","GCF.T1.56","GCF.T1.57","GCF.T1.58","GCF.T1.59","GCF.T1.60","GCF.T1.61","GCF.T1.62","GCF.T1.63","GCF.T1.64","GCF.T1.65","MST.T1.01","MST.T1.02","MST.T1.03","MST.T1.04","MST.T1.05","MST.T1.06","MST.T1.09","MST.T1.10","MST.T1.11","MST.T1.13","MST.T1.14","MST.T1.16","MST.T1.17","MST.T1.18","MST.T1.19","MST.T1.21","MST.T1.22","MST.T1.23","MST.T1.24","MST.T1.25","MST.T1.26","MST.T1.27","MST.T1.28","MST.T1.29","MST.T1.30","MST.T1.31","MST.T1.32","MST.T1.33","MST.T1.34","MST.T1.35","MST.T1.36","MST.T1.39","MST.T1.40","MST.T1.41","MST.T1.42","MST.T1.44","MST.T1.45","MST.T1.46","MST.T1.47","MST.T1.49","MST.T1.50","MST.T1.51","MST.T1.52","MST.T1.55","MST.T1.56","MST.T1.57","MST.T1.58","MST.T1.59","MST.T1.60","MST.T1.61","MST.T1.62","MST.T1.63","MST.T1.64","MST.T1.65")
mbSet<-UpdateSampleItems(mbSet);
mbSet<-ApplyAbundanceFilter(mbSet, "prevalence", 4, 0.2);
mbSet<-ApplyVarianceFilter(mbSet, "iqr", 0.1);
mbSet<-PerformNormalization(mbSet, "none", "colsum", "none", "true");
mbSet<-PrepareCorrExpValues(mbSet, "Source", "Genus", "dbgr", "reingold-tilford", "all", "0.05")
mbSet<-PerformNetworkCorrelation(mbSet,"Genus", "secom_p1", "expr",100, 0.05, 0.3, "mean", "cor_net_1.json")
mbSet<-PlotBoxDataCorr(mbSet, "box_plot_corr_0","Coprococcus", "png", "");
mbSet<-PrepareCorrExpValues(mbSet, "Source", "Genus", "dbgr", "reingold-tilford", "all", "0.05")
mbSet<-PerformNetworkCorrelation(mbSet,"Genus", "secom_p1", "expr",200, 0.05, 0.3, "mean", "cor_net_2.json")
mbSet<-PlotBoxDataCorr(mbSet, "box_plot_corr_1","Coprococcus", "png", "")
