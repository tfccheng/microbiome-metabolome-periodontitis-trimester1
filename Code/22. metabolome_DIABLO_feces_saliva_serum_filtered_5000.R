## -------------------------------------------------------------------------------------------------------------------
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
df_metab_feces_X.log.pareto <- readRDS("Data prepare/df_metab_feces_X.log.pareto.rds")


# load metabolite identification
df_metab_feces_identification <- readRDS("Data filtering/df_metab_feces_identification.rds")


dim(df_metab_feces_X.log.pareto)
metab.meta <- read_excel("metabolome sample meta.xlsx")

#load microbiome


df.micro.feces.g <- readRDS("Data filtering/df.micro.feces.g.rds")
df.micro.plaque.g <- readRDS("Data filtering/df.micro.plaque.g.rds")


#load meta
meta.metab.micro.combine.feces <- readRDS("Data filtering/meta.metab.micro.combine.feces.rds")


####clinical data preparation####
## clinical continuous data
df.clinical.num.impute <- readRDS("Data prepare/df.clinical.num.impute.rds")
df.clinical.num.impute.2 <- readRDS("Data prepare/df.clinical.num.impute.2.rds")
set.seed(1234)



####DIABLO genus level####



#####subgingival plaque microbiome and clinical parameters####
DIABLO.data.plaque.2.g = list(microb = df.micro.plaque.g, 
                            clinic = df.clinical.num.impute.2)
DIABLO.Y <- meta.metab.micro.combine.feces$Disease #identical as feces
lapply(DIABLO.data.plaque.2.g, dim) 
summary(DIABLO.Y)

#initial analysis
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(24, 24)

# generate PLS models
pls.DIABLO.plaque.2.g <- spls(DIABLO.data.plaque.2[["microb"]], DIABLO.data.plaque.2.g[["clinic"]], 
                            keepX = list.keepX, keepY = list.keepY)


# plot features of first PLS
plotVar(pls.DIABLO.plaque.2.g, cutoff = 0.5, title = "microbiome vs clinical", 
        legend = c("microbiome", "clinical"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('cyan4', 'goldenrod2'))



cor(pls.DIABLO.plaque.2.g$variates$X, pls.DIABLO.plaque.2.g$variates$Y) 

#check rownames whether are identical
identical(row.names(df.micro.plaque.g), row.names(df.clinical.num.impute))

design.plaque.2.g = matrix(0.1, ncol = length(DIABLO.data.plaque.2.g), nrow = length(DIABLO.data.plaque.2.g), 
                         dimnames = list(names(DIABLO.data.plaque.2.g), names(DIABLO.data.plaque.2.g)))
diag(design.plaque.2.g) = 0 # set diagonal to 0s

design.plaque.2.g

# form basic DIABLO model
basic.diablo.model.plaque.2.g = block.splsda(X = DIABLO.data.plaque.2.g, Y = DIABLO.Y, ncomp = 5, design = design.plaque.2.g)
perf.diablo.plaque.2.g = perf(basic.diablo.model.plaque.2.g, validation = 'Mfold', 
                            folds = 10, nrepeat = 50, cpus = 12, progressBar = TRUE) 
plot(perf.diablo.plaque.2.g)

ncomp.diablo.plaque.2.g = perf.diablo.plaque.2.g$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp.diablo.plaque.2.g = 2
perf.diablo.plaque.2.g$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX.plaque.2.g = list (microb = c(seq(10, 147, 10)), 
                            clinic = c(seq(4, 18, 2), seq(20,24,3)))
# run the feature selection tuning
tune.diablo.plaque.2.g = tune.block.splsda(X = DIABLO.data.plaque.2.g, Y = DIABLO.Y, ncomp = ncomp.diablo.plaque.2.g, 
                                         test.keepX = test.keepX.plaque.2.g, design = design.plaque.2.g,
                                         validation = 'Mfold', folds = 10, nrepeat = 50, 
                                         BPPARAM = BiocParallel::SnowParam(workers = parallel::detectCores()-2, progressbar = TRUE),
                                         progressBar = TRUE,
                                         dist = "centroids.dist")

list.keepX.plaque.2.g = tune.diablo.plaque.2.g$choice.keepX # set the optimal values of features to retain
tune.diablo.plaque.2.g$choice.ncomp$ncomp #1
list.keepX.plaque.2.g
final.diablo.model.plaque.2.g = block.splsda(X = DIABLO.data.plaque.2.g, Y = DIABLO.Y, ncomp = ncomp.diablo.plaque.2.g, 
                                           keepX = list.keepX.plaque.2.g, design = design.plaque.2.g)

final.diablo.model.plaque.2.g$design
selectVar(final.diablo.model.plaque.2.g, block = 'microb', comp = 1)$microb$name

plotDiablo(final.diablo.model.plaque.2.g, ncomp = 1, col.per.group = c("#00468B", "#AD002A"))

plotIndiv.diablo.plaque2.g <-plotIndiv(final.diablo.model.plaque.2.g, ind.names = FALSE, legend = TRUE,  
                                     col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                     title = 'DIABLO subgingival plaque plots')
pl.plotIndiv.diablo.plaque2.g <- plotIndiv.diablo.plaque2.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotArrow.diablo.plaque2.g <-plotArrow(final.diablo.model.plaque.2.g, ind.names = FALSE, legend = TRUE, 
                                     col.per.group = c("#00468B", "#AD002A"), legend.title = "",
                                     title = 'DIABLO subgingival plaque')
pl.plotArrow.diablo.plaque2.g <- plotArrow.diablo.plaque2.g +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =0.5),
        aspect.ratio = 1)
par(pty="s")
plotVar(final.diablo.model.plaque.2.g, var.names = FALSE, cutoff = 0.5, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 15), cex = c(2,2), 
        col = c('cyan4', 'goldenrod2'))
dev.off()
df.circosPlot.plaque.2.g<-circosPlot(final.diablo.model.plaque.2.g, cutoff = 0.3, line = TRUE, comp = 1,
                                   color.blocks= c('cyan4', 'goldenrod2'),
                                   color.Y = c("#00468B", "#AD002A"),
                                   color.cor = c("#b2182b","#2166ac"), size.variables = 0.75, size.labels = 1.5)
saveRDS(df.circosPlot.plaque.2.g, "DIABLO/plaque_microb.genus.clinic/df.circosPlot.plaque.2.g.rds")

network.plaque.2.g <- network(final.diablo.model.plaque.2.g, blocks = c(1,2), comp = list(microb = 1, clinic = 1),
                            color.node = c('cyan4', 'goldenrod2'), cutoff = 0.3)
write.graph(network.plaque.2.g$gR, file = "DIABLO/plaque_microb.genus.clinic/network.plaque.2.g.gml", format = "gml")

plotLoadings(final.diablo.model.plaque.2.g, comp = 1, contrib = 'max', method = 'median', title = "plaque", ndisplay = 100)

plotLoadings.final.diablo.model.plaque.2.g.microb <- plotLoadings(final.diablo.model.plaque.2.g, comp = 1, block = 1, contrib = 'max', method = 'median', title = "plaque")
df.pld.diablo.plaque.2.g.microb <- plotLoadings.final.diablo.model.plaque.2.g.microb %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
pld.diablo.plaque.2.g.microb <- ggplot(data=df.pld.diablo.plaque.2.g.microb, aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Subgingival plaque\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.diablo.model.plaque.2.g.clinic <- plotLoadings(final.diablo.model.plaque.2.g, comp = 1, block = 2, contrib = 'max', method = 'median', title = "plaque")
df.pld.diablo.plaque.2.g.clinic <- plotLoadings.final.diablo.model.plaque.2.g.clinic %>%
  tibble::rownames_to_column("parameter") %>% 
  filter(GroupContrib != "tie")
pld.diablo.plaque.2.g.clinic <- ggplot(data=df.pld.diablo.plaque.2.g.clinic, aes(x=importance, y = reorder(parameter, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Subgingival plaque\nclinical", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.6), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.cimDiablo.plaque.2.g <- cimDiablo(final.diablo.model.plaque.2.g, title = "DIABLO subgingival plaque", transpose = T,
                                   color.Y = c("#00468B", "#AD002A"), comp = 1,
                                   color.blocks = c('cyan4', 'goldenrod2'))
saveRDS(df.cimDiablo.plaque.2.g, "DIABLO/plaque_microb.genus.clinic/df.cimDiablo.plaque.2.g.rds")
auc.plaque.2.g.microb = auroc(final.diablo.model.plaque.2.g, roc.block = "microb",roc.comp = 1)
auc.plaque.2.g.clinic = auroc(final.diablo.model.plaque.2.g, roc.block = "clinic",roc.comp = 1)

auc.plaque.2.g.microb.data <-auc.plaque.2.g.microb$graph.microb$comp1$data[,c("Specificity","Sensitivity")]
auc.plaque.2.g.clinic.data <-auc.plaque.2.g.clinic$graph.clinic$comp1$data[,c("Specificity","Sensitivity")]
auc.plaque.2.g.microb.data[,"ome"] <- "Microbiome"
auc.plaque.2.g.clinic.data[,"ome"] <- "Clinical"
auc.plaque.2.g.data.combine <- rbind(auc.plaque.2.g.microb.data, auc.plaque.2.g.clinic.data)
auc.plaque.2.g.data.combine$ome <- factor(auc.plaque.2.g.data.combine$ome, 
                                        levels = c("Microbiome", "Clinical"))
pl.auc.plaque.2.g.data.combine<-ggplot(data=auc.plaque.2.g.data.combine, aes(x=Specificity, y=Sensitivity, color=ome)) +
  xlab("100 - Specificity (%)") + ylab("Sensitivity (%)") + 
  geom_line(linewidth=1)+
  geom_abline(intercept = 0, lty = 3, lwd = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+
  scale_color_manual(labels = c(paste("Microbiome:\n", format(auc.plaque.2.g.microb$microb$comp1[1], digits = 3), 
                                      " (", format(auc.plaque.2.g.microb$microb$comp1[2], digits = 2), ")", sep = ""),
                                paste("Clinical:\n", format(auc.plaque.2.g.clinic$clinic$comp1[1], digits = 3), 
                                      " (", format(auc.plaque.2.g.clinic$clinic$comp1[2], digits = 2), ")", sep = "")),
                     values = c('cyan4', 'goldenrod2'))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title= element_blank(),
        legend.position = c(0.8, 0.2),
        aspect.ratio = 1,
        axis.line = element_line(colour = "black"))



#####feces microbiome, clinical and feces metabolome####
DIABLO.data.feces.g = list(microb = df.micro.feces.g, 
                            clinic = df.clinical.num.impute.2,
                            metab = df_metab_feces_X.log.pareto)
DIABLO.Y <- meta.metab.micro.combine.feces$Disease
lapply(DIABLO.data.feces.g, dim) 
summary(DIABLO.Y)

#initial analysis
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1.DIABLO.feces.g <- spls(DIABLO.data.feces.g[["microb"]], DIABLO.data.feces.g[["clinic"]], 
                             keepX = list.keepX, keepY = list.keepY) 
pls2.DIABLO.feces.g <- spls(DIABLO.data.feces.g[["microb"]], DIABLO.data.feces.g[["metab"]], 
                             keepX = list.keepX, keepY = list.keepY)
pls3.DIABLO.feces.g <- spls(DIABLO.data.feces.g[["clinic"]], DIABLO.data.feces.g[["metab"]], 
                             keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1.DIABLO.feces.g, cutoff = 0.4, title = "(a) microbiome vs clinical", 
        legend = c("microbiome", "clinical"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2.DIABLO.feces.g, cutoff = 0.5, title = "(b) microbiome vs metabolome", 
        legend = c("microbiome", "metabolome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3.DIABLO.feces.g, cutoff = 0.5, title = "(c) clinical vs metabolome", 
        legend = c("clinical", "metabolome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

cor(pls1.DIABLO.feces.g$variates$X, pls1.DIABLO.feces.g$variates$Y) 
cor(pls2.DIABLO.feces.g$variates$X, pls2.DIABLO.feces.g$variates$Y)
cor(pls3.DIABLO.feces.g$variates$X, pls3.DIABLO.feces.g$variates$Y) 
#check rownames whether are identical
identical(row.names(df.clinical.num.impute), row.names(df_metab_feces_X.log.pareto))
identical(row.names(df.clinical.num.impute), row.names(df.micro.feces.g))

design.feces.g = matrix(0.1, ncol = length(DIABLO.data.feces.g), nrow = length(DIABLO.data.feces.g), 
                         dimnames = list(names(DIABLO.data.feces.g), names(DIABLO.data.feces.g)))
diag(design.feces.g) = 0 # set diagonal to 0s


design.feces.g

# form basic DIABLO model
basic.diablo.model.feces.g = block.splsda(X = DIABLO.data.feces.g, Y = DIABLO.Y, ncomp = 5, design = design.feces.g)
perf.diablo.feces.g = perf(basic.diablo.model.feces.g, validation = 'Mfold', 
                            folds = 10, nrepeat = 50, cpus = 12, progressBar = TRUE) 
plot(perf.diablo.feces.g)

ncomp.diablo.feces.g = perf.diablo.feces.g$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
perf.diablo.feces.g$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX.feces.g = list (microb = c(seq(10, 150, 10)), 
                            clinic = c(seq(4, 18, 2), seq(20,24,3)),
                            metab = c(seq(150,250,50), seq(300, 1000, 100)))

# run the feature selection tuning
tune.diablo.feces.g = tune.block.splsda(X = DIABLO.data.feces.g, Y = DIABLO.Y, ncomp = ncomp.diablo.feces.g, 
                                         test.keepX = test.keepX.feces.g, design = design.feces.g,
                                         validation = 'Mfold', folds = 10, nrepeat = 50, 
                                         BPPARAM = BiocParallel::SnowParam(workers = parallel::detectCores()-2, progressbar = TRUE),
                                         progressBar = TRUE,
                                         dist = "centroids.dist")   
saveRDS(tune.diablo.feces, "tune.diablo.feces.rds")
# tune.diablo.feces <- readRDS("tune.diablo.feces.rds")


list.keepX.feces.g = tune.diablo.feces.g$choice.keepX # set the optimal values of features to retain
list.keepX.feces.g
tune.diablo.feces.g$choice.ncomp$ncomp #1
final.diablo.model.feces.g = block.splsda(X = DIABLO.data.feces.g, Y = DIABLO.Y, ncomp = ncomp.diablo.feces.g, 
                                           keepX = list.keepX.feces.g, design = design.feces.g)

final.diablo.model.feces.g$design
selectVar(final.diablo.model.feces.g, block = 'microb', comp = 1)$microb$name

plotDiablo(final.diablo.model.feces.g, ncomp = 1, col.per.group = c("#00468B", "#AD002A"))

plotIndiv.diablo.feces.g <- plotIndiv(final.diablo.model.feces.g, ind.names = FALSE, legend = TRUE,
                                    col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                    title = 'DIABLO feces plots')
pl.plotIndiv.diablo.feces.g <- plotIndiv.diablo.feces.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotArrow.diablo.feces.g <-plotArrow(final.diablo.model.feces.g, ind.names = FALSE, legend = TRUE, 
                                   col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                   title = 'DIABLO feces')
pl.plotArrow.diablo.feces.g <- plotArrow.diablo.feces.g +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =0.5),
        aspect.ratio = 1)

par(pty="s")
plotVar(final.diablo.model.feces.g, var.names = FALSE, cutoff = 0.5, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 15, 17), cex = c(2,2,2), 
        col = c('cyan4', 'goldenrod2', 'indianred3'))
dev.off()

df.circosPlot.feces.g <- circosPlot(final.diablo.model.feces.g, cutoff = 0.3, line = TRUE, comp = 1,
                                  color.blocks= c('cyan4', 'goldenrod2', 'indianred3'),
                                  color.Y = c("#00468B", "#AD002A"),
                                  color.cor = c("#b2182b","#2166ac"), size.labels = 1.5)
saveRDS(df.circosPlot.feces.g, "DIABLO/feces_microb.genus.clinic.metab/df.circosPlot.feces.g.rds")
network.feces.g <- network(final.diablo.model.feces.g, blocks = c(1,2,3), comp = list(microb =1, clinic =1, metab =1),
                         color.node = c('cyan4', 'goldenrod2', 'indianred3'), cutoff = 0.3)
write.graph(network.feces.g$gR, file = "DIABLO/feces_microb.genus.clinic.metab/network.feces.g.gml", format = "gml")

plotLoadings(final.diablo.model.feces.g, comp = 1, contrib = 'max', method = 'median', title = "feces", ndisplay = 40)

plotLoadings.final.diablo.model.feces.g.microb <- plotLoadings(final.diablo.model.feces.g, comp = 1, block = 1, contrib = 'max', method = 'median', title = "feces")
df.pld.diablo.feces.g.microb <- plotLoadings.final.diablo.model.feces.g.microb %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
pld.diablo.feces.g.microb <- ggplot(data=df.pld.diablo.feces.g.microb, aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.8, 0.8), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.diablo.model.feces.g.clinic <- plotLoadings(final.diablo.model.feces.g, comp = 1, block = 2, contrib = 'max', method = 'median', title = "feces")
df.pld.diablo.feces.g.clinic <- plotLoadings.final.diablo.model.feces.g.clinic %>%
  tibble::rownames_to_column("parameter") %>%
  filter(GroupContrib != "tie")
pld.diablo.feces.g.clinic <- ggplot(data=df.pld.diablo.feces.g.clinic, aes(x=importance, y = reorder(parameter, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nclinical", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.diablo.model.feces.g.metab <- plotLoadings(final.diablo.model.feces.g, comp = 1, block = 3, contrib = 'max', method = 'median', title = "feces")

df.pld.diablo.feces.g.metab <- plotLoadings.final.diablo.model.feces.g.metab %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_feces_identification, by = "ID") 

df.pld.diablo.feces.g.metab$identification_pld <- ifelse(is.na(df.pld.diablo.feces.g.metab$identification), 
                                                       df.pld.diablo.feces.g.metab$ID, 
                                                       df.pld.diablo.feces.g.metab$identification)

df.pld.diablo.feces.g.metab.identified <- df.pld.diablo.feces.g.metab %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.g.metab.identified$identification_id <- str_trunc(df.pld.diablo.feces.g.metab.identified$identification_id, 100)

pld.diablo.feces.g.metab.identified <- ggplot(data=df.pld.diablo.feces.g.metab.identified, 
                                            aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.diablo.feces.g.metab.MS2 <- df.pld.diablo.feces.g.metab %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.g.metab.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                         df.pld.diablo.feces.g.metab.MS2$identification_MS2)
pld.diablo.feces.g.metab.MS2 <- ggplot(data=df.pld.diablo.feces.g.metab.MS2 %>% filter(abs(importance)>0.05),
                                     aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                         fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.diablo.feces.g.metab.single.id <- df.pld.diablo.feces.g.metab %>%
  filter(!is.na(identification)) %>% #filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.diablo.feces.g.metab.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                     df.pld.diablo.feces.g.metab.single.id$identification_single_id)
pld.diablo.feces.g.metab.single.id <- ggplot(data=df.pld.diablo.feces.g.metab.single.id, 
                                           aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                               fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.cimDiablo.feces.g <- cimDiablo(final.diablo.model.feces.g, title = "DIABLO feces", transpose = T,
                                  color.Y = c("#00468B", "#AD002A"), comp = 1,
                                  color.blocks = c('cyan4', 'goldenrod2', 'indianred3'))
saveRDS(df.cimDiablo.feces.g, "DIABLO/feces_microb.genus.clinic.metab/df.cimDiablo.feces.g.rds")
auc.feces.g.microb = auroc(final.diablo.model.feces.g, roc.block = "microb", roc.comp = 1)
auc.feces.g.clinic = auroc(final.diablo.model.feces.g, roc.block = "clinic", roc.comp = 1)
auc.feces.g.metab = auroc(final.diablo.model.feces.g, roc.block = "metab", roc.comp = 1)
auc.feces.g.microb.data <-auc.feces.g.microb$graph.microb$comp1$data[,c("Specificity","Sensitivity")]
auc.feces.g.clinic.data <-auc.feces.g.clinic$graph.clinic$comp1$data[,c("Specificity","Sensitivity")]
auc.feces.g.metab.data <-auc.feces.g.metab$graph.metab$comp1$data[,c("Specificity","Sensitivity")]
auc.feces.g.microb.data[,"ome"] <- "Microbiome"
auc.feces.g.clinic.data[,"ome"] <- "Clinical"
auc.feces.g.metab.data[,"ome"] <- "Metabolome"
auc.feces.g.data.combine <- rbind(auc.feces.g.microb.data, auc.feces.g.clinic.data, auc.feces.g.metab.data)
auc.feces.g.data.combine$ome <- factor(auc.feces.g.data.combine$ome, 
                                       levels = c("Microbiome", "Clinical", "Metabolome"))
pl.auc.feces.g.data.combine<-ggplot(data=auc.feces.g.data.combine, aes(x=Specificity, y=Sensitivity, color=ome)) +
  xlab("100 - Specificity (%)") + ylab("Sensitivity (%)") + 
  geom_line(linewidth=1)+
  geom_abline(intercept = 0, lty = 3, lwd = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+ 
  scale_y_continuous(expand = c(0, 0), limits = c(-5,105), n.breaks = 6)+
  scale_color_manual(labels = c(paste("Microbiome:\n", format(auc.feces.g.microb$microb$comp1[1], digits = 3), 
                                      " (", format(auc.feces.g.microb$microb$comp1[2], digits = 2), ")", sep = ""),
                                paste("Clinical:\n", format(auc.feces.g.microb$clinic$comp1[1], digits = 3), 
                                      " (", format(auc.feces.g.microb$clinic$comp1[2], digits = 2), ")", sep = ""),
                                paste("Metabolome:\n", format(auc.feces.g.microb$metab$comp1[1], digits = 3), 
                                      " (", format(auc.feces.g.microb$metab$comp1[2], digits = 2), ")", sep = "")),
                     values = c('cyan4', 'goldenrod2', 'indianred3'))+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title= element_blank(),
        legend.position = c(0.8, 0.2),
        aspect.ratio = 1,
        axis.line = element_line(colour = "black"))


####edit nodes of network####
#import metabolite pathway data
idms2.kegg.feces <- read_excel("Metab identification/feces/idms2.kegg_metabolism_pathway.xlsx")
idms2.kegg.saliva <- read_excel("Metab identification/saliva/idms2.kegg_metabolism_pathway.xlsx")
idms2.kegg.serum <- read_excel("Metab identification/serum/idms2.kegg_metabolism_pathway.xlsx")

#####feces microb clinic metab genus####
network.feces.g.gml.nodes <- read.csv("DIABLO/network.feces.g.gml default node.csv") %>%
  mutate(Genus = label) %>%
  mutate(ID = label) %>%
  mutate(parameter = label) %>%
  dplyr::left_join(data.frame(rowData(altExp(tse.micro.subsamples[["MST"]], "Genus"))), by = "Genus") %>%
  dplyr::left_join(df.pld.diablo.feces.g.metab.single.id, by = "ID") %>%
  dplyr::left_join(df.pld.diablo.feces.g.microb %>% dplyr::rename(Genus = ASV) %>% 
                     dplyr::select(Genus, GroupContrib), by = "Genus") %>%
  dplyr::left_join(df.pld.diablo.feces.clinic %>% 
                     dplyr::select(parameter, GroupContrib), by = "parameter")

#add pathway info
network.feces.g.gml.nodes$pathway <- NA
for (i in 1:dim(network.feces.g.gml.nodes)[1]){
  tmp <- grepl(network.feces.g.gml.nodes$label[i], idms2.kegg.feces$Feature)
  network.feces.g.gml.nodes$pathway[i] <-  idms2.kegg.feces[tmp,]$Level2[1]
}
write.csv(network.feces.g.gml.nodes,
          "DIABLO/network.feces.g.gml.nodes.csv", row.names = F)


#####plaque genus micro clinic####
network.plaque.2.g.gml.nodes <- read.csv("DIABLO/network.plaque.2.g.gml default node.csv") %>%
  mutate(Genus = label) %>%
  mutate(parameter = label) %>%
  dplyr::left_join(data.frame(rowData(altExp(tse.micro.subsamples[["SBP"]], "Genus"))), by = "Genus") %>%
  dplyr::left_join(df.pld.diablo.plaque.2.g.microb %>% dplyr::rename(Genus = ASV) %>% 
                     dplyr::select(Genus, GroupContrib), by = "Genus") %>%
  dplyr::left_join(df.pld.diablo.plaque.2.g.clinic %>% 
                     dplyr::select(parameter, GroupContrib), by = "parameter")
write.csv(network.plaque.2.g.gml.nodes,
          "DIABLO/network.plaque.2.g.gml.nodes.csv", row.names = F)


