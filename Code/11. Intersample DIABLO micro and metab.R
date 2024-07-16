library(mixOmics) # import the mixOmics library
library(dplyr)
library(readxl)
library(mia)
library(stringr)
library(ggplot2)
library(igraph)

setwd("")

####metabolome data preparation####
# load IQR filtered metabolome
df_metab_feces_X.log.pareto <- readRDS("")
df_metab_saliva_X.log.pareto <- readRDS("")
df_metab_serum_X.log.pareto <- readRDS("")
# load metabolite identification
df_metab_feces_identification <- readRDS("")
df_metab_saliva_identification <- readRDS("")
df_metab_serum_identification <- readRDS("")

#load microbiome


df.micro.saliva.g <- readRDS("")
df.micro.feces.g <- readRDS("")
df.micro.plaque.g <- readRDS("")



#load meta
meta.metab.micro.combine.saliva <- readRDS("")


#####metablome in saliva, serum, feces####
DIABLO.data.metab = list(saliva = df_metab_saliva_X.log.pareto, 
                          serum = df_metab_serum_X.log.pareto,
                          feces = df_metab_feces_X.log.pareto)
DIABLO.Y <- meta.metab.micro.combine.saliva$Disease
lapply(DIABLO.data.metab, dim) 
summary(DIABLO.Y)

#initial analysis
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1.DIABLO.metab <- spls(DIABLO.data.metab[["saliva"]], DIABLO.data.metab[["serum"]], 
                           keepX = list.keepX, keepY = list.keepY) 
pls2.DIABLO.metab <- spls(DIABLO.data.metab[["saliva"]], DIABLO.data.metab[["feces"]], 
                           keepX = list.keepX, keepY = list.keepY)
pls3.DIABLO.metab <- spls(DIABLO.data.metab[["serum"]], DIABLO.data.metab[["feces"]], 
                           keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1.DIABLO.metab, cutoff = 0.5, title = "(a) saliva vs serum", 
        legend = c("saliva", "serum"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('cyan4', 'goldenrod2'))

# plot features of second PLS
plotVar(pls2.DIABLO.metab, cutoff = 0.5, title = "(b) saliva vs feces", 
        legend = c("saliva", "feces"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('cyan4', 'indianred3'))

# plot features of third PLS
plotVar(pls3.DIABLO.metab, cutoff = 0.5, title = "(c) serum vs feces", 
        legend = c("serum", "feces"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('goldenrod2', 'indianred3'))

cor(pls1.DIABLO.metab$variates$X, pls1.DIABLO.metab$variates$Y) 
cor(pls2.DIABLO.metab$variates$X, pls2.DIABLO.metab$variates$Y)
cor(pls3.DIABLO.metab$variates$X, pls3.DIABLO.metab$variates$Y) 
#check rownames whether are identical
identical(row.names(df_metab_saliva_X.log.pareto), row.names(df_metab_serum_X.log.pareto))
identical(row.names(df_metab_saliva_X.log.pareto), row.names(df_metab_feces_X.log.pareto))

design.metab = matrix(0.1, ncol = length(DIABLO.data.metab), nrow = length(DIABLO.data.metab), 
                       dimnames = list(names(DIABLO.data.metab), names(DIABLO.data.metab)))
diag(design.metab) = 0 # set diagonal to 0s

# design.metab[1,3] = 0.5
# design.metab[3,1] = 0.5
design.metab

# form basic DIABLO model
basic.diablo.model.metab = block.splsda(X = DIABLO.data.metab, Y = DIABLO.Y, ncomp = 5, design = design.metab)
perf.diablo.metab = perf(basic.diablo.model.metab, validation = 'Mfold', 
                          folds = 10, nrepeat = 50, cpus = 12, progressBar = TRUE) 
plot(perf.diablo.metab)

ncomp.diablo.metab = perf.diablo.metab$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
perf.diablo.metab$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX.metab = list (saliva = c(seq(10, 100, 20), seq(150,250,50), seq(300, 500, 100)), 
                         serum = c(seq(10, 100, 20), seq(150,250,50), seq(300, 500, 100)),
                         feces = c(seq(10, 100, 20), seq(150,250,50), seq(300, 500, 100)))

# run the feature selection tuning
tune.diablo.metab = tune.block.splsda(X = DIABLO.data.metab, Y = DIABLO.Y, ncomp = ncomp.diablo.metab, 
                                       test.keepX = test.keepX.metab, design = design.metab,
                                       validation = 'Mfold', folds = 10, nrepeat = 50, 
                                       BPPARAM = BiocParallel::SnowParam(workers = parallel::detectCores()-2, progressbar = TRUE),
                                       progressBar = TRUE,
                                       dist = "centroids.dist")



list.keepX.metab = tune.diablo.metab$choice.keepX # set the optimal values of features to retain
list.keepX.metab
tune.diablo.metab$choice.ncomp$ncomp #1
final.diablo.model.metab = block.splsda(X = DIABLO.data.metab, Y = DIABLO.Y, ncomp = ncomp.diablo.metab, 
                                         keepX = list.keepX.metab, design = design.metab)

final.diablo.model.metab$design
selectVar(final.diablo.model.metab, block = 'feces', comp = 1)$feces$name

plotDiablo(final.diablo.model.metab, ncomp = 1, col.per.group = c("#00468B", "#AD002A"))

plotIndiv.diablo.metab <- plotIndiv(final.diablo.model.metab, ind.names = FALSE, legend = TRUE,
                                     col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                     title = 'DIABLO metab plots')
pl.plotIndiv.diablo.metab <- plotIndiv.diablo.metab$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotArrow.diablo.metab <-plotArrow(final.diablo.model.metab, ind.names = FALSE, legend = TRUE, 
                                    col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                    title = 'DIABLO metab')
pl.plotArrow.diablo.metab <- plotArrow.diablo.metab +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =0.5),
        aspect.ratio = 1)

par(pty="s")
plotVar(final.diablo.model.metab, var.names = FALSE, cutoff = 0.5, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 15, 17), cex = c(2,2,2), 
        col = c('cyan4', 'goldenrod2', 'indianred3'))
dev.off()

df.circosPlot.metab <- circosPlot(final.diablo.model.metab, cutoff = 0.5, line = TRUE, comp = 1,
                                   color.blocks= c('cyan4', 'goldenrod2', 'indianred3'),
                                   color.Y = c("#00468B", "#AD002A"),
                                   color.cor = c("#b2182b","#2166ac"), size.labels = 1.5)
saveRDS(df.circosPlot.metab, "DIABLO/intersample.metab/df.circosPlot.metab.rds")

network.metab <- network(final.diablo.model.metab, blocks = c(1,2,3),
                          color.node = c('cyan4', 'goldenrod2', 'indianred3'), cutoff = 0.3)
write.graph(network.metab$gR, file = "DIABLO/intersample.metab/network.metab.gml", format = "gml")


plotLoadings(final.diablo.model.metab, comp = 1, contrib = 'max', method = 'median', title = "metab", ndisplay = 40)


#saliva
plotLoadings.final.diablo.model.saliva.metab <- plotLoadings(final.diablo.model.metab, comp = 1, block = 1, contrib = 'max', method = 'median', title = "saliva")

df.pld.diablo.saliva.metab <- plotLoadings.final.diablo.model.saliva.metab %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_saliva_identification, by = "ID") 

df.pld.diablo.saliva.metab$identification_pld <- ifelse(is.na(df.pld.diablo.saliva.metab$identification), 
                                                        df.pld.diablo.saliva.metab$ID, 
                                                        df.pld.diablo.saliva.metab$identification)

df.pld.diablo.saliva.metab.identified <- df.pld.diablo.saliva.metab %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.saliva.metab.identified$identification_id <- str_trunc(df.pld.diablo.saliva.metab.identified$identification_id, 100)

pld.diablo.saliva.metab.identified <- ggplot(data=df.pld.diablo.saliva.metab.identified, 
                                             aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                 fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.1), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.diablo.saliva.metab.MS2 <- df.pld.diablo.saliva.metab %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.saliva.metab.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                          df.pld.diablo.saliva.metab.MS2$identification_MS2)
pld.diablo.saliva.metab.MS2 <- ggplot(data=df.pld.diablo.saliva.metab.MS2, 
                                      aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                          fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.diablo.saliva.metab.single.id <- df.pld.diablo.saliva.metab %>%
  filter(!is.na(identification)) %>% #filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.saliva.metab.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                      df.pld.diablo.saliva.metab.single.id$identification_single_id)
pld.diablo.saliva.metab.single.id <- ggplot(data=df.pld.diablo.saliva.metab.single.id, 
                                            aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#serum
plotLoadings.final.diablo.model.serum.metab <- plotLoadings(final.diablo.model.metab, comp = 1, block = 2, contrib = 'max', method = 'median', title = "serum")

df.pld.diablo.serum.metab <- plotLoadings.final.diablo.model.serum.metab %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_serum_identification, by = "ID") 

df.pld.diablo.serum.metab$identification_pld <- ifelse(is.na(df.pld.diablo.serum.metab$identification), 
                                                        df.pld.diablo.serum.metab$ID, 
                                                        df.pld.diablo.serum.metab$identification)

df.pld.diablo.serum.metab.identified <- df.pld.diablo.serum.metab %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.serum.metab.identified$identification_id <- str_trunc(df.pld.diablo.serum.metab.identified$identification_id, 100)

pld.diablo.serum.metab.identified <- ggplot(data=df.pld.diablo.serum.metab.identified, 
                                             aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                 fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.1), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.diablo.serum.metab.MS2 <- df.pld.diablo.serum.metab %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.serum.metab.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                          df.pld.diablo.serum.metab.MS2$identification_MS2)
pld.diablo.serum.metab.MS2 <- ggplot(data=df.pld.diablo.serum.metab.MS2, 
                                      aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                          fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.6, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.diablo.serum.metab.single.id <- df.pld.diablo.serum.metab %>%
  filter(!is.na(identification)) %>% #filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.serum.metab.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                      df.pld.diablo.serum.metab.single.id$identification_single_id)
pld.diablo.serum.metab.single.id <- ggplot(data=df.pld.diablo.serum.metab.single.id, 
                                            aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.6, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#feces
plotLoadings.final.diablo.model.feces.metab <- plotLoadings(final.diablo.model.metab, comp = 1, block = 3, contrib = 'max', method = 'median', title = "feces")

df.pld.diablo.feces.metab <- plotLoadings.final.diablo.model.feces.metab %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_feces_identification, by = "ID") 

df.pld.diablo.feces.metab$identification_pld <- ifelse(is.na(df.pld.diablo.feces.metab$identification), 
                                                        df.pld.diablo.feces.metab$ID, 
                                                        df.pld.diablo.feces.metab$identification)

df.pld.diablo.feces.metab.identified <- df.pld.diablo.feces.metab %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.metab.identified$identification_id <- str_trunc(df.pld.diablo.feces.metab.identified$identification_id, 100)

pld.diablo.feces.metab.identified <- ggplot(data=df.pld.diablo.feces.metab.identified, 
                                             aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                 fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.1), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.diablo.feces.metab.MS2 <- df.pld.diablo.feces.metab %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.metab.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                          df.pld.diablo.feces.metab.MS2$identification_MS2)
pld.diablo.feces.metab.MS2 <- ggplot(data=df.pld.diablo.feces.metab.MS2, 
                                      aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                          fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.5, 0.5), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.diablo.feces.metab.single.id <- df.pld.diablo.feces.metab %>%
  filter(!is.na(identification)) %>% #filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.metab.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                      df.pld.diablo.feces.metab.single.id$identification_single_id)
pld.diablo.feces.metab.single.id <- ggplot(data=df.pld.diablo.feces.metab.single.id, 
                                            aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.5, 0.5), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.cimDiablo.metab <- cimDiablo(final.diablo.model.metab, title = "DIABLO metab", transpose = T,
                                 color.Y = c("#00468B", "#AD002A"), 
                                 color.blocks = c('cyan4', 'goldenrod2', 'indianred3'))
saveRDS(df.cimDiablo.metab, "DIABLO/intersample.metab/df.cimDiablo.metab.rds")
auc.metab.saliva = auroc(final.diablo.model.metab, roc.block = "saliva", roc.comp = 1)
auc.metab.serum = auroc(final.diablo.model.metab, roc.block = "serum", roc.comp = 1)
auc.metab.feces = auroc(final.diablo.model.metab, roc.block = "feces", roc.comp = 1)
auc.metab.saliva.data <-auc.metab.saliva$graph.saliva$comp1$data[,c("Specificity","Sensitivity")]
auc.metab.serum.data <-auc.metab.serum$graph.serum$comp1$data[,c("Specificity","Sensitivity")]
auc.metab.feces.data <-auc.metab.feces$graph.feces$comp1$data[,c("Specificity","Sensitivity")]
auc.metab.saliva.data[,"ome"] <- "Saliva"
auc.metab.serum.data[,"ome"] <- "Serum"
auc.metab.feces.data[,"ome"] <- "Feces"
auc.metab.data.combine <- rbind(auc.metab.saliva.data, auc.metab.serum.data, auc.metab.feces.data)
auc.metab.data.combine$ome <- factor(auc.metab.data.combine$ome, 
                                      levels = c("Saliva", "Serum", "Feces"))
pl.auc.metab.data.combine<-ggplot(data=auc.metab.data.combine, aes(x=Specificity, y=Sensitivity, color=ome)) +
  xlab("100 - Specificity (%)") + ylab("Sensitivity (%)") +
  geom_line(linewidth=1)+
  geom_abline(intercept = 0, lty = 3, lwd = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+
  scale_color_manual(labels = c(paste("Saliva:\n", format(auc.metab.saliva$saliva$comp1[1], digits = 3), 
                                      " (", format(auc.metab.saliva$saliva$comp1[2], digits = 2), ")", sep = ""),
                                paste("Serum:\n", format(auc.metab.serum$serum$comp1[1], digits = 3), 
                                      " (", format(auc.metab.serum$serum$comp1[2], digits = 2), ")", sep = ""),
                                paste("Feces:\n", format(auc.metab.feces$feces$comp1[1], digits = 3), 
                                      " (", format(auc.metab.feces$feces$comp1[2], digits = 2), ")", sep = "")),
                     values = c('cyan4', 'goldenrod2', 'indianred3'))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title= element_blank(),
        legend.position = c(0.8, 0.2),
        aspect.ratio = 1,
        axis.line = element_line(colour = "black"))







#####microbiome genus in plaque, saliva, feces####
DIABLO.data.microb.g = list(saliva = df.micro.saliva.g, 
                          plaque = df.micro.plaque.g,
                          feces = df.micro.feces.g)
DIABLO.Y <- meta.metab.micro.combine.saliva$Disease
lapply(DIABLO.data.microb.g, dim) 
summary(DIABLO.Y)

#initial analysis
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1.DIABLO.microb.g <- spls(DIABLO.data.microb.g[["saliva"]], DIABLO.data.microb.g[["plaque"]], 
                           keepX = list.keepX, keepY = list.keepY) 
pls2.DIABLO.microb.g <- spls(DIABLO.data.microb.g[["saliva"]], DIABLO.data.microb.g[["feces"]], 
                           keepX = list.keepX, keepY = list.keepY)
pls3.DIABLO.microb.g <- spls(DIABLO.data.microb.g[["plaque"]], DIABLO.data.microb.g[["feces"]], 
                           keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1.DIABLO.microb.g, cutoff = 0.5, title = "(a) saliva vs plaque", 
        legend = c("saliva", "plaque"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('cyan4', 'goldenrod2'))

# plot features of second PLS
plotVar(pls2.DIABLO.microb.g, cutoff = 0.5, title = "(b) saliva vs feces", 
        legend = c("saliva", "feces"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('cyan4', 'indianred3'))

# plot features of third PLS
plotVar(pls3.DIABLO.microb.g, cutoff = 0.5, title = "(c) plaque vs feces", 
        legend = c("plaque", "feces"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('goldenrod2', 'indianred3'))

cor(pls1.DIABLO.microb.g$variates$X, pls1.DIABLO.microb.g$variates$Y) 
cor(pls2.DIABLO.microb.g$variates$X, pls2.DIABLO.microb.g$variates$Y)
cor(pls3.DIABLO.microb.g$variates$X, pls3.DIABLO.microb.g$variates$Y) 
#check rownames whether are identical
identical(row.names(df.micro.saliva), row.names(df.micro.plaque))
identical(row.names(df.micro.saliva), row.names(df.micro.feces))

design.microb.g = matrix(0.1, ncol = length(DIABLO.data.microb.g), nrow = length(DIABLO.data.microb.g), 
                       dimnames = list(names(DIABLO.data.microb.g), names(DIABLO.data.microb.g)))
diag(design.microb.g) = 0 # set diagonal to 0s

# design.microb.g[1,3] = 0.5
# design.microb.g[3,1] = 0.5
design.microb.g

# form basic DIABLO model
basic.diablo.model.microb.g = block.splsda(X = DIABLO.data.microb.g, Y = DIABLO.Y, ncomp = 5, design = design.microb.g)
perf.diablo.microb.g = perf(basic.diablo.model.microb.g, validation = 'Mfold', 
                          folds = 10, nrepeat = 50, cpus = 12, progressBar = TRUE) 
plot(perf.diablo.microb.g)

ncomp.diablo.microb.g = perf.diablo.microb.g$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
perf.diablo.microb.g$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX.microb.g = list (saliva = c(seq(10, 120, 10)), 
                          plaque = c(seq(10, 120, 10)),
                          feces = c(seq(10, 120, 10)))

# run the feature selection tuning
tune.diablo.microb.g = tune.block.splsda(X = DIABLO.data.microb.g, Y = DIABLO.Y, ncomp = ncomp.diablo.microb.g, 
                                       test.keepX = test.keepX.microb.g, design = design.microb.g,
                                       validation = 'Mfold', folds = 10, nrepeat = 50, 
                                       BPPARAM = BiocParallel::SnowParam(workers = parallel::detectCores()-2, progressbar = TRUE),
                                       progressBar = TRUE,
                                       dist = "centroids.dist")
list.keepX.microb.g = tune.diablo.microb.g$choice.keepX # set the optimal values of features to retain
list.keepX.microb.g
tune.diablo.microb.g$choice.ncomp$ncomp #2


final.diablo.model.microb.g = block.splsda(X = DIABLO.data.microb.g, Y = DIABLO.Y, ncomp = ncomp.diablo.microb.g, 
                                        keepX = list.keepX.microb.g, design = design.microb.g)

final.diablo.model.microb.g$design
selectVar(final.diablo.model.microb.g, block = 'feces', comp = 1)$feces$name

plotDiablo(final.diablo.model.microb.g, ncomp = 1, col.per.group = c("#00468B", "#AD002A"))

plotIndiv.diablo.microb.g <- plotIndiv(final.diablo.model.microb.g, ind.names = FALSE, legend = TRUE,
                                    col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                    title = 'DIABLO microb plots')
pl.plotIndiv.diablo.microb.g <- plotIndiv.diablo.microb.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotArrow.diablo.microb.g <-plotArrow(final.diablo.model.microb.g, ind.names = FALSE, legend = TRUE, 
                                   col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                   title = 'DIABLO microb')
pl.plotArrow.diablo.microb.g <- plotArrow.diablo.microb.g +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =0.5),
        aspect.ratio = 1)

par(pty="s")
plotVar(final.diablo.model.microb.g, var.names = FALSE, cutoff = 0.5, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 15, 17), cex = c(2,2,2), 
        col = c('cyan4', 'goldenrod2', 'indianred3'))
dev.off()

df.circosPlot.microb.g <- circosPlot(final.diablo.model.microb.g, cutoff = 0.3, line = TRUE, comp = 1,
                                  color.blocks= c('cyan4', 'goldenrod2', 'indianred3'),
                                  color.Y = c("#00468B", "#AD002A"),
                                  color.cor = c("#b2182b","#2166ac"), size.labels = 1.5)
saveRDS(df.circosPlot.microb.g, "DIABLO/intersample.microb.g/df.circosPlot.microb.g.rds")

network.microb.g <- network(final.diablo.model.microb.g, blocks = c(1,2,3),
                         color.node = c('cyan4', 'goldenrod2', 'indianred3'), cutoff = 0.3)
write.graph(network.microb.g$gR, file = "DIABLO/intersample.microb.g/network.microb.g.gml", format = "gml")


plotLoadings(final.diablo.model.microb.g, comp = 1, contrib = 'max', method = 'median', title = "microb", ndisplay = 40)
plotLoadings.final.diablo.model.microb.g.saliva <- plotLoadings(final.diablo.model.microb.g, comp = 1, block = 1, contrib = 'max', method = 'median', title = "saliva")
df.pld.diablo.microb.g.saliva <- plotLoadings.final.diablo.model.microb.g.saliva %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
pld.diablo.microb.g.saliva <- ggplot(data=df.pld.diablo.microb.g.saliva, aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.diablo.model.microb.g.plaque <- plotLoadings(final.diablo.model.microb.g, comp = 1, block = 2, contrib = 'max', method = 'median', title = "plaque")
df.pld.diablo.microb.g.plaque <- plotLoadings.final.diablo.model.microb.g.plaque %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
pld.diablo.microb.g.plaque <- ggplot(data=df.pld.diablo.microb.g.plaque, aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Plaque\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.diablo.model.microb.g.feces <- plotLoadings(final.diablo.model.microb.g, comp = 1, block = 3, contrib = 'max', method = 'median', title = "feces")
df.pld.diablo.microb.g.feces <- plotLoadings.final.diablo.model.microb.g.feces %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
pld.diablo.microb.g.feces <- ggplot(data=df.pld.diablo.microb.g.feces, aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.5, 0.5), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.cimDiablo.microb.g <- cimDiablo(final.diablo.model.microb.g, title = "DIABLO microb", transpose = T,
                                color.Y = c("#00468B", "#AD002A"), 
                                color.blocks = c('cyan4', 'goldenrod2', 'indianred3'))
saveRDS(df.cimDiablo.microb.g, "DIABLO/intersample.microb.g/df.cimDiablo.microb.g.rds")
auc.microb.g.saliva = auroc(final.diablo.model.microb.g, roc.block = "saliva", roc.comp = 1)
auc.microb.g.plaque = auroc(final.diablo.model.microb.g, roc.block = "plaque", roc.comp = 1)
auc.microb.g.feces = auroc(final.diablo.model.microb.g, roc.block = "feces", roc.comp = 1)
auc.microb.g.saliva.data <-auc.microb.g.saliva$graph.saliva$comp1$data[,c("Specificity","Sensitivity")]
auc.microb.g.plaque.data <-auc.microb.g.plaque$graph.plaque$comp1$data[,c("Specificity","Sensitivity")]
auc.microb.g.feces.data <-auc.microb.g.feces$graph.feces$comp1$data[,c("Specificity","Sensitivity")]
auc.microb.g.saliva.data[,"ome"] <- "Saliva"
auc.microb.g.plaque.data[,"ome"] <- "Plaque"
auc.microb.g.feces.data[,"ome"] <- "Feces"
auc.microb.g.data.combine <- rbind(auc.microb.g.saliva.data, auc.microb.g.plaque.data, auc.microb.g.feces.data)
auc.microb.g.data.combine$ome <- factor(auc.microb.g.data.combine$ome, 
                                     levels = c("Saliva", "Plaque", "Feces"))
pl.auc.microb.g.data.combine<-ggplot(data=auc.microb.g.data.combine, aes(x=Specificity, y=Sensitivity, color=ome)) +
  xlab("100 - Specificity (%)") + ylab("Sensitivity (%)") +
  geom_line(linewidth=1)+
  geom_abline(intercept = 0, lty = 3, lwd = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+
  scale_color_manual(labels = c(paste("Saliva:\n", format(auc.microb.g.saliva$saliva$comp1[1], digits = 3), 
                                      " (", format(auc.microb.g.saliva$saliva$comp1[2], digits = 2), ")", sep = ""),
                                paste("Subgingival plaque:\n", format(auc.microb.g.plaque$plaque$comp1[1], digits = 3), 
                                      " (", format(auc.microb.g.plaque$plaque$comp1[2], digits = 2), ")", sep = ""),
                                paste("Feces:\n", format(auc.microb.g.feces$feces$comp1[1], digits = 3), 
                                      " (", format(auc.microb.g.feces$feces$comp1[2], digits = 2), ")", sep = "")),
                     values = c('cyan4', 'goldenrod2', 'indianred3'))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title= element_blank(),
        legend.position = c(0.8, 0.2),
        aspect.ratio = 1,
        axis.line = element_line(colour = "black"))


