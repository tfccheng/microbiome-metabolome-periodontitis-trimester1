library(microbiome)
library(phyloseq)
library(biomformat)
library(dplyr)

set.seed(123)
setwd("")

#load prefiltered data
pseq.silva.minimal.filter <- readRDS("")

pseq.silva.minimal.filter %>%
  refseq() %>%
  Biostrings::writeXStringSet("PICRUSt2/asv.minimal.filter.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
asv_tab_for_picrust2.minimal.filter <- as.data.frame(otu_table(pseq.silva.minimal.filter))
write.table(asv_tab_for_picrust2.minimal.filter, 
            "PICRUSt2/asv_tab_for_picrust2.minimal.filter.txt", sep = "\t", quote = F, col.names = NA)



#filter data of features with 5% sample with 1 count
rmn_feat <- ncol(otu_table(pseq.silva.minimal.filter))
smpl.perc <- 0.05 #percentage of samples
cnt <- 1 #minimal counts of a feature
minLen <- smpl.perc*rmn_feat
kept.inx <- apply(otu_table(pseq.silva.minimal.filted), 
                  MARGIN = 1, function(x) {sum(x >= cnt) >= minLen})
pseq.silva.prune.0.05.prev <- prune_taxa(kept.inx, pseq.silva.minimal.filted) 

pseq.silva.prune.0.05.prev %>%
  refseq() %>%
  Biostrings::writeXStringSet("PICRUSt2/asv.prune.0.05.prev.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

asv_tab_for_picrust2.prune.0.05.prev <- as.data.frame(otu_table(pseq.silva.prune.0.05.prev))
write.table(asv_tab_for_picrust2.prune.0.05.prev, 
            "PICRUSt2/asv_tab_for_picrust2.prune.0.05.prev.txt", sep = "\t", quote = F, col.names = NA)



