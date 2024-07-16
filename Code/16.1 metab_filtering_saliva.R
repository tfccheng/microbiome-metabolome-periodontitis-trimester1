library(MetaboAnalystR)
setwd("Data filtering/saliva/")
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "saliva_combine-norm-metaboAnalystInput_PDNPD.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "iqr", 40, "T", 20, F)
mSet<-PreparePrenormData(mSet)
mSet<-GetGroupNames(mSet, "")
feature.nm.vec <- c("")
smpl.nm.vec <- c("")
grp.nm.vec <- c("NPD","PD")
mSet<-UpdateData(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, T, 0.05, FALSE, FALSE, "raw", FALSE) #wilcox, p<0.05
mSet<-PlotTT(mSet, "wilcox.saliva.p0.05", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, FALSE, "raw", FALSE) #t-test, p<0.05
mSet<-PlotTT(mSet, "tt.saliva.p0.05", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, T, 0.2, FALSE, FALSE, "fdr", FALSE)  #wilcox, FDR<0.2
mSet<-PlotTT(mSet, "wilcox.saliva.fdr0.2", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.2, FALSE, FALSE, "fdr", FALSE)  #t-test, FDR<0.2
mSet<-PlotTT(mSet, "tt.saliva.fdr0.2", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.2, FALSE, FALSE, "fdr", FALSE)

mSet<-Ttests.Anal(mSet, F, 1, FALSE, FALSE, "raw", FALSE) #t-test, keep all

#volcano
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.05, FALSE, "raw")  #t-test, p<0.05
mSet<-PlotVolcano(mSet, "volcano.tt.saliva.p0.05.FC2.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.2, FALSE, "fdr")  #t-test, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.tt.saliva.fdr0.2.FC2.",1, 0, "pdf", 72, width=NA)

mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, T, 0.05, FALSE, "raw")  #wilcox, p<0.05
mSet<-PlotVolcano(mSet, "volcano.wilcox.saliva.p0.05.FC2.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, T, 0.2, FALSE, "fdr")  #wilcox, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.wilcox.saliva.fdr0.2.FC2.",1, 0, "pdf", 72, width=NA)

#volcano FC =1.5
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.05, FALSE, "raw")  #t-test, p<0.05
mSet<-PlotVolcano(mSet, "volcano.tt.serum.p0.05.FC1.5.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.2, FALSE, "fdr")  #t-test, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.tt.serum.fdr0.2.FC1.",1, 0, "pdf", 72, width=NA)

mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, T, 0.05, FALSE, "raw")  #wilcox, p<0.05
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.p0.05.FC1.5.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, T, 0.2, FALSE, "fdr")  #wilcox, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.fdr0.2.FC1.5.",1, 0, "pdf", 72, width=NA)

#SAM analysis parametric
mSet<-SAM.Anal(mSet, "d.stat", FALSE, FALSE, 0.59, "sam_saliva_para.")
mSet<-PlotSAM.Cmpd(mSet, "sam_saliva_para.", "pdf", 72, width=NA)
mSet<-PlotSAM.FDR(mSet, "sam_saliva_para.FDR.", "pdf", 72, width=NA)

mSet<-SAM.Anal(mSet, "wilc.stat", FALSE, FALSE, 0.22, "sam_saliva_nonpara.")
mSet<-PlotSAM.Cmpd(mSet, "sam_saliva_nonpara.", "pdf", 72, width=NA)
mSet<-PlotSAM.FDR(mSet, "sam_saliva_nonpara.FDR.", "pdf", 72, width=NA)

