library(microbiome)
library(mia)
library(phyloseq)
library(knitr)
library(dplyr)
library(readxl)
library(ggplot2)

setwd("")


##### Step 1 incorporate data into phyloseq Class object with minimal filter (singleton) #####

set.seed(123)
asv_tab <- read.csv("ASVs_counts.txt", sep = "\t", row.names = 1) %>% as.matrix()
tax_tab.silva <- read.csv("ASV_silva_species_combined_genus.txt", sep = "\t", row.names = 1) %>% as.matrix()
meta_tab <- read.csv("meta new.csv", row.names = 1)
phy_tre <- read_tree("T1.tre")
ref_seq <- Biostrings::readDNAStringSet("ASVs.fa", format="fasta")
pseq.silva <- phyloseq(otu_table(asv_tab, taxa_are_rows = T),
                       tax_table(tax_tab.silva),
                       sample_data(meta_tab),
                       phy_tree(phy_tre),
                       refseq(ref_seq))
saveRDS(pseq.silva, "pseq.silva.rds")
tse.raw <- makeTreeSummarizedExperimentFromPhyloseq(pseq.silva)
scater::perCellQCMetrics(tse.raw)
tse.raw <- scater::addPerCellQC(tse.raw)
p1 <- ggplot(colData(tse.raw)) +
  geom_histogram(aes(x = sum), color = "black", 
                 fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  theme_bw() +
  # Removes the grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # Adds y-axis
        axis.line = element_line(colour = "black")) 
df <- as.data.frame(colData(tse.raw)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df, aes(y = index, x = sum/1e6, color = Sample)) +
  geom_point() +  
  labs(x = "Library size (million reads)", 
       y = "Sample index") +  
  scale_color_manual(values = c(SBP="#00a087", MSA="#0f77c1", MST ="#e64b35")) +
  theme_bw() +
  # Removes the grid
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # Adds y-axis
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1)


rare_color <- c(rep("#00a087",54), rep("#0f77c1",54), rep("#e64b35",54))
vegan::rarecurve(t(assay(tse.raw,"counts")[]), 3496, ylab = "Species richness", label = F, col = rare_color, cex = 0.6) 

pseq.silva.minimal.filter <- filter_taxa(pseq.silva, function (x) {sum(x > 0) > 1}, prune=TRUE) 
saveRDS(pseq.silva.minimal.filter, "pseq.silva.minimal.filter.rds")
pseq.silva.minimal.filter <- readRDS("pseq.silva.minimal.filter.rds")


pseq.silva.minimal.filter.f <- microbiome::add_besthit(pseq.silva.minimal.filter)
taxa_names(pseq.silva.minimal.filter.f) <- 
  gsub(":.*\\.", ":", taxa_names(pseq.silva.minimal.filter.f)) #clean names
saveRDS(pseq.silva.minimal.filter.f, "pseq.silva.minimal.filter.f.rds")

pseq.silva.minimal.filter.f <- readRDS("pseq.silva.minimal.filter.f.rds")



##### Step 2 convert phyloseq object to  TreeSummarizedExperiment object#####
tse.silva <- makeTreeSummarizedExperimentFromPhyloseq(pseq.silva.minimal.filter.f)
saveRDS(tse.silva, "tse.silva.rds")
clinicaldataID <- rownames(colData(tse.silva))
colData(tse.silva) <- cbind(colData(tse.silva), clinicaldataID)
colData(tse.silva)$clinicaldataID <- gsub("[0-9A-Za-z]+.T1.", "T1-", colData(tse.silva)$clinicaldataID)
#subset by Samples
tse.silva.SBP <- tse.silva[ , tse.silva$Sample %in% c("SBP")]
tse.silva.MSA <- tse.silva[ , tse.silva$Sample %in% c("MSA")]
tse.silva.MST <- tse.silva[ , tse.silva$Sample %in% c("MST")]

saveRDS(tse.silva.SBP, "tse.silva.SBP.rds")
saveRDS(tse.silva.MSA, "tse.silva.MSA.rds")
saveRDS(tse.silva.MST, "tse.silva,MST.rds")

# clinical data
df.clinical <- read_excel("T1 data 55 subjects with 16s clinical data.xlsx", na = "NA") 
df.clinical.final <- df.clinical[!(df.clinical$clinicaldataID =="T1-53"),] %>% # remove subject T1-53, not be sequenced
  rename(DiseaseGroup = Disease) # rename column Disease
saveRDS(df.clinical.final, "T1clinical.rds")

#cbind colData with clinical data
#check the order of ID before cbind
tse.silva.wclinical <- tse.silva
colData(tse.silva.wclinical) <- cbind(colData(tse.silva.wclinical), df.clinical.final)
saveRDS(tse.silva.wclinical, "tse.silva.wclinical.rds")


tse.silva.SBP.wclinical <- tse.silva.SBP
tse.silva.MSA.wclinical <- tse.silva.MSA
tse.silva.MST.wclinical <- tse.silva.MST

colData(tse.silva.SBP.wclinical) <- cbind(colData(tse.silva.SBP.wclinical), df.clinical.final)
colData(tse.silva.MSA.wclinical) <- cbind(colData(tse.silva.MSA.wclinical), df.clinical.final)
colData(tse.silva.MST.wclinical) <- cbind(colData(tse.silva.MST.wclinical), df.clinical.final)



saveRDS(tse.silva.SBP.wclinical, "tse.silva.SBP.wclinical.rds")
saveRDS(tse.silva.MSA.wclinical, "tse.silva.MSA.wclinical.rds")
saveRDS(tse.silva.MST.wclinical, "tse.silva,MST.wclinical.rds")

tse <- tse.silva.wclinical

##clean clinical data
#gestationalWeek
tse$gestationalWeek[c(2, 56, 110)] <- NA
#MenstrualDays use average for range data
tse$MenstrualDays[c(28, 82, 136)] <- 6
tse$MenstrualDays[c(43, 97, 151)] <- 6.5
#MenstrualCycle use average for range data
tse$MenstrualCycle[c(52, 106, 160)] <- 32.5

# transform these as numeric
numeric.names <- c("gestationalWeek", "MenstrualDays", "MenstrualCycle",
                   "parity","abortionNo")
colData(tse)[, numeric.names] <- lapply(colData(tse)[, numeric.names], as.numeric)

colnames(colData(tse))[14] <- "PD4"
colnames(colData(tse))[22] <- "AAP_EFP_PD"
colnames(colData(tse))[23] <- "AAP_EFP_AL"
#filter no use data
colData(tse) <- cbind(colData(tse)[,c("Source", "Sample", "Disease", "ID", "clinicaldataID", "teethno", "FMPS", "BOP", 
                                      "PD4", "CAL12", "CDC", "AAP_EFP_PD", "AAP_EFP_AL", "FDIgraph",
                                      "FDIscore_total", "FDIscore_classification", "extent", "age",
                                      "education", "income", "fertilization", "gestationalWeek", "pregnantBefore",
                                      "pastGPA_G", "pastGPA_P", "pastGPA_A", 
                                      "menarche", "regularMenstruation", "MenstrualDays", "MenstrualCycle")],
                      colData(tse)[,56:109])  #from perioTreatment  
#transform these as categorical
category.names <- c("Source", "Sample", "Disease","CDC", "AAP_EFP_PD", "AAP_EFP_AL", "FDIgraph",
                    "FDIscore_total", "FDIscore_classification", "extent", "education", "income", "fertilization", "pregnantBefore",
                    "pastGPA_G", "pastGPA_P", "pastGPA_A", "regularMenstruation", "perioTreatment", 
                    "regularDentalVisit","perioFamilyHistory", "teethbrushing", "BOB", "dentalFloss", "otherdentalcleaning",
                    "smoking", "passiveSmoking", "drinking", "diet", "saltedFood", "smokedFood", "vegeFruit",              
                    "cornsBeans", "highCarboFatFood", "probiotics", "yogurt", "exercise", "exerciseTime", "thyroidfunction")
colData(tse)[, category.names] <-lapply(colData(tse)[, category.names], factor)


df.caprotectin <- read.csv("calprotectin-finalest_w_ID.csv")
a <- data.frame(colData(tse))
a <- a %>%
  dplyr::left_join(df.caprotectin[,c(1,3)], by = "ID")
colData(tse) <-cbind(colData(tse), a$calprotectin)

saveRDS(tse, "tse.silva.wclean.clinical.calprotectin.rds")
