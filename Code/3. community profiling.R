library(microbiome)
library(phyloseq)
library(knitr)
library(dplyr)
library(vegan)
library(ggplot2)
library(mia)
library(miaViz)
library(scater)

set.seed(123)
setwd("")


pseq.silva.minimal.filted <- readRDS("pseq.silva.minimal.filter.rds")
pseq.gtdb.minimal.filted <- readRDS("pseq.gtdb.minimal.filter.rds")

## Use microbiome package
# data transformation
pseq.silva.comptrans <- microbiome::transform(pseq.silva.minimal.filted,
                                              "compositional")
sample_data(pseq.silva.comptrans)$Sample <- factor(sample_data(pseq.silva.comptrans)$Sample, levels = c("SBP", "MSA", "MST"))
sample_data(pseq.silva.comptrans)$Sample <- recode(sample_data(pseq.silva.comptrans)$Sample, 
                                                   SBP = "Subgingival plaque",
                                                   MSA = "Saliva",
                                                   MST = "Feces")

#pseq.silva.comptrans.phylum <- aggregate_taxa(pseq.silva.comptrans,level = "Phylum")
pseq.silva.comptrans.phylum <- aggregate_rare(pseq.silva.comptrans,
                                             level = "Phylum",
                                             detection = 1/100,
                                             prevalence = 1/100)
pseq.silva.comptrans.genus <- aggregate_rare(pseq.silva.comptrans,
                                       level = "Genus",
                                       detection = 1/100,
                                       prevalence = 10/100)
pseq.silva.comptrans.family <- aggregate_rare(pseq.silva.comptrans,
                                             level = "Family",
                                             detection = 1/100,
                                             prevalence = 10/100)
# Limit the analysis on core taxa and specific sample group
library(hrbrthemes)
library(RColorBrewer)
#library(pals)
mycolor <- colorRampPalette(brewer.pal(12, "Paired"))(30)
#mycolor <- colorRampPalette(brewer.pal(8, "Set1"))(32)
mycolor1 <- colorRampPalette(brewer.pal(12, "Paired"))(13)
p.rel.abund.all.phylum <- plot_composition(pseq.silva.comptrans.phylum,
                      taxonomic.level = "Phylum",
                      group_by = "Sample", sample.sort = "Bacteroidota",
                      otu.sort = "abundance") +
  scale_fill_manual("Phylum", values = mycolor1) +
  labs(y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom")
p.rel.abund.all.phylum

p.rel.abund.all.phylum.avebysource <- plot_composition(pseq.silva.comptrans.phylum,
                                                       otu.sort = "abundance",
                                                       group_by = "Sample",
                                                       average_by = "Source") +
  scale_fill_manual("Phylum", values = mycolor1) +
  labs(y = "Relative abundance") +
  theme_ipsum(grid= F,
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.title.y = element_text(hjust = 0.5))
p.rel.abund.all.phylum.avebysource
write.csv(p.rel.abund.all.phylum.avebysource$data, "Figures/1 Composition/p.rel.abund.all.phylum.avebysource.csv")

p.rel.abund.all.genus <- plot_composition(pseq.silva.comptrans.genus,
                                          taxonomic.level = "Genus",
                                          group_by = "Sample",
                                          otu.sort = "abundance") +
  scale_fill_manual("Genus", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(face = "italic"))
p.rel.abund.all.genus

p.rel.abund.all.genus.avebysource <- plot_composition(pseq.silva.comptrans.genus,
                                                      otu.sort = "abundance",
                                                      group_by = "Sample",
                                                      average_by = "Source") +
  scale_fill_manual("Genus", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.title.y = element_text(hjust = 0.5),
        legend.text = element_text(face = "italic"))

p.rel.abund.all.genus.avebysource
write.csv(p.rel.abund.all.genus.avebysource$data, "Figures/1 Composition/p.rel.abund.all.genus.avebysource.csv")

p.rel.abund.all.family <- plot_composition(pseq.silva.comptrans.family,
                                          taxonomic.level = "Family",
                                          group_by = "Sample",
                                          otu.sort = "abundance") +
  scale_fill_manual("Family", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom")
p.rel.abund.all.family

p.rel.abund.all.family.avebysource <- plot_composition(pseq.silva.comptrans.family,
                                                      otu.sort = "abundance",
                                                      group_by = "Sample",
                                                      average_by = "Source") +
  scale_fill_manual("Family", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F, 
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.title.y = element_text(hjust = 0.5))
p.rel.abund.all.family.avebysource
write.csv(p.rel.abund.all.family.avebysource$data, "Figures/1 Composition/p.rel.abund.all.family.avebysource.csv")

