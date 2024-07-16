library(dada2)
library(ggplot2)
library(dplyr)

set.seed(100)
#change according Windows or Linux
setwd("")
seq_folder <- "./RawData" #change to the directory containing the fastq raw files
analysis_path <- "./R"


list.files(seq_folder)
file_folders <-list.files(seq_folder)
my_16S_folder_path <- vector(length(file_folders), mode="character")
fnFs <- vector(length(file_folders), mode="character")
fnRs <- vector(length(file_folders), mode="character")

for (i in 1:length(file_folders)) {
  my_16S_folder_path[i] <- paste(seq_folder, file_folders[i], sep="/")
  fnFs[i] <- grep(list.files(my_16S_folder_path[i], pattern="[0-9A-Za-z]_1.f(q|astq).gz", recursive = TRUE, full.names = TRUE), 
                  pattern='raw_', invert=TRUE, value=TRUE)
  fnRs[i] <- grep(list.files(my_16S_folder_path[i], pattern="[0-9A-Za-z]_2.f(q|astq).gz", recursive = TRUE, full.names = TRUE), 
                  pattern='raw_', invert=TRUE, value=TRUE)
}

sample.names <- file_folders
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

 #use MicrobiomeAnalystR data2_utilities.R for dada2
#run data2_utilities_revised.R first

#Windows
setParametersRes <- setParameters(file_compressed = TRUE, 
              OS_is_windows = TRUE)
#Linux
setParametersRes <- setParameters(file_compressed = TRUE, 
                                  OS_is_windows = FALSE)

.plotQualityProfileLoop(f_F = fnFs, #forward reads file
                        f_R = fnRs,#reverse reads file
                        sn = sample.names, # sample name
                        plot_format = "pdf", # pdf, tiff, ....
                        fd = analysis_path)

#DON'T use the function seAanityCheck, because raw files are not in the same folder
#use the following script instead
seqSanityCheckRes <- list(fnFs = fnFs,
                          fnRs = fnRs,
                          sn = sample.names)


processRawSeqRes <- processRawSeq(setParametersRes = setParametersRes, # results from setParameters
              seqSanityCheckRes = seqSanityCheckRes, # results from seqSanityCheck
              reads_trim_length_F_R = '', #for V3-V4 no need for this setting #retained reads length for forward and reverse,
              #depending on the quailty of reads, the quailty graph can be obtained by seqSanityCheck results;
              plot_format = "pdf")

saveRDS(processRawSeqRes, "processRawSeqRes.rds")
constructSeqTabRes <- constructSeqTab(setParametersRes = setParametersRes, # results from setParameters
                                      processRawSeqRes = processRawSeqRes)

saveRDS(constructSeqTabRes, "constructSeqTabRes.rds")

#check the dataset URL first in data2_utilities.R
options(timeout = 1000) # set timeout limit for download.file
assignTaxRes <- assignTax(constructSeqTabRes = constructSeqTabRes, #results from constructSeqTab
          ref_db = "silva")
saveRDS(assignTaxRes, "assignTaxRes.rds")
#assignTaxRes_gtdb <- assignTax(constructSeqTabRes = constructSeqTabRes, #results from constructSeqTab
#                          ref_db = "gtdb")



# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(constructSeqTabRes$seqtab.nochim)
asv_headers <- vector(dim(constructSeqTabRes$seqtab.nochim)[2], mode="character")
asv_headers_not_fasta <- vector(dim(constructSeqTabRes$seqtab.nochim)[2], mode="character")
for (i in 1:dim(constructSeqTabRes$seqtab.nochim)[2]) {
  asv_headers_not_fasta[i] <- paste("ASV", i, sep="_")
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "sequence_table/ASVs.fa")

# count table:
asv_tab <- t(constructSeqTabRes$seqtab.nochim)
row.names(asv_tab) <- asv_headers_not_fasta
write.table(cbind.data.frame("#NAME" = row.names(asv_tab),
                             asv_tab),
            file = file.path("sequence_table/ASVs_counts.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)




#taxonomy
asv_tax_silva_spp <- as.data.frame(read.csv("tax/ASV_silva_species.txt", sep = "\t"))
asv_tax_silva_spp <- cbind.data.frame("#TAXONOMY" = asv_headers_not_fasta,
                                      asv_tax_silva_spp) %>%
  dplyr::select(-2)

# rename Prevotella_7 and Prevotella_9 as Prevotella
library(stringi)
asv_tax_silva_spp$Genus1 <- stri_replace_all_regex(asv_tax_silva_spp$Genus,
                                                   pattern = c("_7", "_9", " group"), 
                                                   replacement = "", 
                                                   vectorize = FALSE)
asv_tax_silva_spp <- asv_tax_silva_spp %>%
  mutate(Species = if_else(!is.na(Species) & !is.na(Genus1), 
                           paste(Genus1, Species, sep = "_"), 
                           Species)) %>%
  dplyr::select(-Genus1)


write.table(asv_tax_silva_spp,
            file = "tax/ASV_silva_species_combined_genus.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")



#export biom file

library(biomformat)
asv_tab_for_biom <- as.data.frame(asv_tab)
biom16S <- make_biom(asv_tab_for_biom)
biom_file <- biom16S
outfile = "sequence_table/OTU_table.biom"
write_biom(biom_file, outfile)

# for phylogenetic tree
library(phangorn)
library(DECIPHER)
seqs <- getSequences(constructSeqTabRes$seqtab.nochim)
names(seqs) <- asv_headers_not_fasta # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=TRUE, processors = NULL)
saveRDS(alignment, "alignment.rds") 
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
write.phyDat(phangAlign, file="phangAlign.fasta", format="fasta")
saveRDS(phangAlign, "phangAlign.rds")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
write.tree(treeNJ, file="treeNJ.tre")

#use FASTTREE on https://usegalaxy.org/ to generate ML tree



