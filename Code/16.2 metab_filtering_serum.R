library(MetaboAnalystR)
setwd("Data filtering/serum/")
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "serum_combine-norm-metaboAnalystInput_PDNPD.csv", "colu", "disc");
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
mSet<-PlotTT(mSet, "wilcox.serum.p0.05", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, FALSE, "raw", FALSE) #t-test, p<0.05
mSet<-PlotTT(mSet, "tt.serum.p0.05", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, T, 0.2, FALSE, FALSE, "fdr", FALSE)  #wilcox, FDR<0.2
mSet<-PlotTT(mSet, "wilcox.serum.fdr0.2", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.2, FALSE, FALSE, "fdr", FALSE)  #t-test, FDR<0.2
mSet<-PlotTT(mSet, "tt.serum.fdr0.2", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.2, FALSE, FALSE, "fdr", FALSE)

mSet<-Ttests.Anal(mSet, F, 1, FALSE, FALSE, "raw", FALSE) #t-test, keep all

#volcano FC =2
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.05, FALSE, "raw")  #t-test, p<0.05
mSet<-PlotVolcano(mSet, "volcano.tt.serum.p0.05.FC2.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.2, FALSE, "fdr")  #t-test, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.tt.serum.fdr0.2.FC2.",1, 0, "pdf", 72, width=NA)

mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, T, 0.05, FALSE, "raw")  #wilcox, p<0.05
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.p0.05.FC2.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, T, 0.2, FALSE, "fdr")  #wilcox, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.fdr0.2.FC2.",1, 0, "pdf", 72, width=NA)
#volcano FC =1.5
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.05, FALSE, "raw")  #t-test, p<0.05
mSet<-PlotVolcano(mSet, "volcano.tt.serum.p0.05.FC1.5.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, F, 0.2, FALSE, "fdr")  #t-test, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.tt.serum.fdr0.2.FC1.",1, 0, "pdf", 72, width=NA)

mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, T, 0.05, FALSE, "raw")  #wilcox, p<0.05
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.p0.05.FC1.5.",1, 0, "pdf", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 1.5, 0, T, 0.2, FALSE, "fdr")  #wilcox, fdr<0.2
mSet<-PlotVolcano(mSet, "volcano.wilcox.serum.fdr0.2.FC1.5.",1, 0, "pdf", 72, width=NA)

