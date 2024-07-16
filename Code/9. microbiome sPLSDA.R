library(mixOmics) # import the mixOmics library
library(dplyr)
library(ggpubr)
library(ggrepel)

setwd("")
#load microbiome
df.micro.saliva <- readRDS("")
df.micro.feces <- readRDS("")
df.micro.plaque <- readRDS("")

df.micro.saliva.g <- readRDS("")
df.micro.feces.g <- readRDS("")
df.micro.plaque.g <- readRDS("")

df.micro.saliva.s <- readRDS("")
df.micro.feces.s <- readRDS("")
df.micro.plaque.s <- readRDS("")

dim(df.micro.saliva)
meta <- readRDS("")
Y.sub <- meta$Disease
Y.sub.color <- case_when(Y.sub == "NPD" ~ "#00468B", 
                         Y.sub == "PD" ~ "#AD002A")
taxa.micro <- read.csv("Data prepare/taxa.micro.csv") %>%
  dplyr::rename(ASV = X)
taxa.micro.g <- read.csv("Data prepare/taxa.micro.genus.csv") %>%
  select(-X)
taxa.micro.sp <- read.csv("Data prepare/taxa.micro.species.csv") %>%
  select(-X)

set.seed(1234)

###################################################################
#######################     SPLSDA     ############################
###################################################################
####ASV####
##### subgingival plaque ####
pca.SBP = pca(df.micro.plaque, ncomp = 10, center = T, scale = T) 
plot(pca.SBP)
plotInd.pca.SBP <- plotIndiv(pca.SBP, group = Y.sub, ind.names = FALSE, # plot the samples projected
          legend = FALSE, title = 'Subgingival plaque')
biplot(pca.SBP, group = Y.sub, 
       legend.title = '')
splsda.SBP <- splsda(df.micro.plaque, Y.sub, ncomp = 10)

plotInd.splsda.SBP <- plotIndiv(splsda.SBP, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = FALSE, title = 'Subgingival plaque')

background.SBP = background.predict(splsda.SBP, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotInd.splsda.SBP.background <- plotIndiv(splsda.SBP, comp = 1:2,
          group = Y.sub, ind.names = FALSE, # colour points by class
          background = background.SBP, # include prediction background for each class
          legend = FALSE, title = "Subgingival plaque")

perf.splsda.SBP <- perf(splsda.SBP, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,# use repeated cross-validation
                                progressBar = T, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.SBP, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.SBP$choice.ncomp


list.keepX.splsda.SBP <- c(seq(10, 100, 10), seq(150,250,50), seq(300, 1175, 100))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.SBP <- tune.splsda(df.micro.plaque, Y.sub, ncomp = 2, 
                                       validation = 'Mfold',
                                       folds = 10, nrepeat = 50, # use repeated cross-validation
                                       dist = 'max.dist', # use max.dist measure
                                       measure = "BER", # use balanced error rate of dist measure
                                       test.keepX = list.keepX.splsda.SBP,
                                       progressBar = T,
                                       cpus = 14) 

plot(tune.splsda.SBP, col = color.jet(2))

tune.splsda.SBP$choice.ncomp$ncomp #1
tune.splsda.SBP$choice.keepX

optimal.ncomp.splsda.SBP <- tune.splsda.SBP$choice.ncomp$ncomp
optimal.ncomp.splsda.SBP <- 2
optimal.keepX.splsda.SBP <- tune.splsda.SBP$choice.keepX[1:optimal.ncomp.splsda.SBP]

final.splsda.SBP <- splsda(df.micro.plaque, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.SBP, 
                                   keepX = optimal.keepX.splsda.SBP)

plotLoadings(final.splsda.SBP, comp = 1, contrib = 'max', method = 'median', title = "SBP", ndisplay = 40)
plotLoadings.final.splsda.SBP <- plotLoadings(final.splsda.SBP, comp = 1, contrib = 'max', method = 'median', title = "SBP")
df.pld.final.splsda.SBP <- plotLoadings.final.splsda.SBP %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie") %>%
  left_join(taxa.micro, by = "ASV") # append taxa info
  
saveRDS(df.pld.final.splsda.SBP, "microbiome sPLSDA/df.pld.final.splsda.SBP.rds")
write.csv(df.pld.final.splsda.SBP, "microbiome sPLSDA/df.pld.final.splsda.SBP.csv")
pld.final.splsda.SBP <- ggplot(data=df.pld.final.splsda.SBP, 
                                aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Subgingival plaque\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.SBP <- plotIndiv(final.splsda.SBP, comp = c(1,2), # plot samples from final model
          group = Y.sub, ind.names = FALSE, # colour by class label
          col.per.group = c("#00468B", "#AD002A"),
          ellipse = TRUE, legend = FALSE, # include 95% confidence ellipse
          title = 'Subgingival plaque')
pl.plotInd.final.splsda.SBP <- plotInd.final.splsda.SBP$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

legend=list(legend = levels(Y.sub), # set of classes
            col = c("#00468B", "#AD002A"), # set of colours
            title = "", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim.SBP <-cim(final.splsda.SBP, row.sideColors = Y.sub.color, title = "Subgingival plaque",
              comp = 1,
              transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.SBP.final <- perf(final.splsda.SBP, 
                          folds = 10, nrepeat = 50, cpus = 14,
                          validation = "Mfold", dist = "max.dist", 
                          progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.SBP.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)


par(mfrow=c(1,1))

plotVar(final.splsda.SBP, comp = c(1,2), cex = 3)

auc.splsda.SBP = auroc(final.splsda.SBP, roc.comp = 1, print = FALSE) 
final.splsda.SBP.selectVar.1 <- selectVar(final.splsda.SBP, comp = 1)$name
final.splsda.SBP.selectVar <- final.splsda.SBP.selectVar.1
df.micro.plaque.select <- df.micro.plaque[,final.splsda.SBP.selectVar.1]

##### saliva ####
pca.MSA = pca(df.micro.saliva, ncomp = 10, center = T, scale = T) 
plot(pca.MSA)
plotInd.pca.MSA <- plotIndiv(pca.MSA, group = Y.sub, ind.names = FALSE, 
          legend = FALSE, title = 'Saliva')


splsda.MSA <- splsda(df.micro.saliva, Y.sub, ncomp = 10)

plotInd.splsda.MSA <- plotIndiv(splsda.MSA, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = FALSE, title = 'Saliva')

background.MSA = background.predict(splsda.MSA, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MSA.background <- plotIndiv(splsda.MSA, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.MSA, 
          legend = FALSE, title = "Saliva")

perf.splsda.MSA<- perf(splsda.MSA, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

plot(perf.splsda.MSA, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MSA$choice.ncomp


list.keepX.splsda.MSA <- c(seq(10, 100, 10), seq(150,250,50), seq(300, 1377, 100))


tune.splsda.MSA <- tune.splsda(df.micro.saliva, Y.sub, ncomp = 4, 
                                        validation = 'Mfold',
                                        folds = 10, nrepeat = 50, 
                                        dist = 'max.dist', 
                                        measure = "BER", 
                                        test.keepX = list.keepX.splsda.MSA,
                                        progressBar = T,
                                        cpus = 14) 

plot(tune.splsda.MSA, col = color.jet(4))

tune.splsda.MSA$choice.ncomp$ncomp #2
tune.splsda.MSA$choice.keepX

optimal.ncomp.splsda.MSA <- tune.splsda.MSA$choice.ncomp$ncomp
#optimal.ncomp.splsda.MSA <- 2
optimal.keepX.splsda.MSA <- tune.splsda.MSA$choice.keepX[1:optimal.ncomp.splsda.MSA]

final.splsda.MSA <- splsda(df.micro.saliva, Y.sub, 
                                    ncomp = optimal.ncomp.splsda.MSA, 
                                    keepX = optimal.keepX.splsda.MSA)

plotLoadings.final.splsda.MSA <- plotLoadings(final.splsda.MSA, comp = 1, contrib = 'max', method = 'median', title = "MSA")
df.pld.final.splsda.MSA <- plotLoadings.final.splsda.MSA %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie") %>%
  left_join(taxa.micro, by = "ASV")

saveRDS(df.pld.final.splsda.MSA, "microbiome sPLSDA/df.pld.final.splsda.MSA.rds")
write.csv(df.pld.final.splsda.MSA, "microbiome sPLSDA/df.pld.final.splsda.MSA.csv")
pld.final.splsda.MSA <- ggplot(data=df.pld.final.splsda.MSA, 
                                aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.3), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plotLoadings.final.splsda.MSA.comp2 <- plotLoadings(final.splsda.MSA, comp = 2, contrib = 'max', method = 'median', title = "MSA")
df.pld.final.splsda.MSA.comp2 <- plotLoadings.final.splsda.MSA.comp2 %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
saveRDS(df.pld.final.splsda.MSA.comp2, "microbiome sPLSDA/df.pld.final.splsda.MSA.comp2.rds")
pld.final.splsda.MSA.comp2 <- ggplot(data=df.pld.final.splsda.MSA.comp2, 
                               aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.3), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.MSA <- plotIndiv(final.splsda.MSA, comp = c(1,2), 
          group = Y.sub, ind.names = FALSE, 
          col.per.group = c("#00468B", "#AD002A"),legend.title = "",
          ellipse = TRUE, legend = FALSE, 
          title = 'Saliva')

pl.plotInd.final.splsda.MSA <- plotInd.final.splsda.MSA$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)


cim.MSA <- cim(final.splsda.MSA, row.sideColors = Y.sub.color, title = "Saliva", 
               comp = 1,
               transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.MSA.final <- perf(final.splsda.MSA, 
                                      folds = 10, nrepeat = 50, cpus = 14,
                                      validation = "Mfold", dist = "max.dist", 
                                      progressBar = TRUE)

par(mfrow=c(1,2))
plot(perf.splsda.MSA.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.MSA.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.MSA, comp = c(1,2), cex = 3)

auc.splsda.MSA = auroc(final.splsda.MSA, roc.comp = 1, print = FALSE) 
final.splsda.MSA.selectVar.1 <- selectVar(final.splsda.MSA, comp = 1)$name
final.splsda.MSA.selectVar.2 <- selectVar(final.splsda.MSA, comp = 2)$name
df.micro.saliva.select <- df.micro.saliva[,unique(c(final.splsda.MSA.selectVar.1,
                                             final.splsda.MSA.selectVar.2))]
##### feces ####

pca.MST = pca(df.micro.feces, ncomp = 10, center = TRUE, scale = T) 
plot(pca.MST)
plotInd.pca.MST <- plotIndiv(pca.MST, group = Y.sub, ind.names = FALSE, 
          legend = FALSE, title = 'Feces')


splsda.MST <- splsda(df.micro.feces, Y.sub, ncomp = 10)

plotInd.splsda.MST <- plotIndiv(splsda.MST, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = FALSE, title = 'Feces')

background.MST = background.predict(splsda.MST, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MST.background <- plotIndiv(splsda.MST, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.MST, 
          legend = FALSE, title = "Feces")

perf.splsda.MST <- perf(splsda.MST, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

plot(perf.splsda.MST, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MST$choice.ncomp


list.keepX.splsda.MST <- c(seq(20, 100, 10), seq(150, 1076, 50))
tune.splsda.MST <- tune.splsda(df.micro.feces, Y.sub, ncomp = 2, 
                                       validation = 'Mfold',
                                       folds = 10, nrepeat = 50, 
                                       dist = 'max.dist', 
                                       measure = "BER", 
                                       test.keepX = list.keepX.splsda.MST,
                                       progressBar = T,
                                       cpus = 14) 

plot(tune.splsda.MST, col = color.jet(2))

tune.splsda.MST$choice.ncomp$ncomp 
tune.splsda.MST$choice.keepX

optimal.ncomp.splsda.MST <- tune.splsda.MST$choice.ncomp$ncomp
#optimal.ncomp.splsda.MST <- 2
optimal.keepX.splsda.MST <- tune.splsda.MST$choice.keepX[1:optimal.ncomp.splsda.MST]

final.splsda.MST <- splsda(df.micro.feces, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.MST, 
                                   keepX = optimal.keepX.splsda.MST)
plotLoadings.final.splsda.MST <- plotLoadings(final.splsda.MST, comp = 1, contrib = 'max', method = 'median', title = "MST")
df.pld.final.splsda.MST <- plotLoadings.final.splsda.MST %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")%>%
  left_join(taxa.micro, by = "ASV")

saveRDS(df.pld.final.splsda.MST, "microbiome sPLSDA/df.pld.final.splsda.MST.rds")
write.csv(df.pld.final.splsda.MST, "microbiome sPLSDA/df.pld.final.splsda.MST.csv")
pld.final.splsda.MST <- ggplot(data=df.pld.final.splsda.MST, 
                                aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.splsda.MST.comp2 <- plotLoadings(final.splsda.MST, comp = 2, contrib = 'max', method = 'median', title = "MST")
df.pld.final.splsda.MST.comp2 <- plotLoadings.final.splsda.MST.comp2 %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
saveRDS(df.pld.final.splsda.MST.comp2, "microbiome sPLSDA/df.pld.final.splsda.MST.comp2.rds")

pld.final.splsda.MST.comp2 <- ggplot(data=df.pld.final.splsda.MST.comp2, 
                               aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.MST <- plotIndiv(final.splsda.MST, comp = c(1,2), 
                                      group = Y.sub, ind.names = FALSE, 
                                      col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                      ellipse = TRUE, legend = FALSE, 
                                      title = 'Feces')
pl.plotInd.final.splsda.MST <- plotInd.final.splsda.MST$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

cim.MST <- cim(final.splsda.MST, row.sideColors = Y.sub.color, title = "Feces",
               comp = 1, 
               transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)
perf.splsda.MST.final <- perf(final.splsda.MST, 
                                folds = 10, nrepeat = 50, cpus = 14,
                                validation = "Mfold", dist = "max.dist", 
                                progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.MST.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.MST.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.MST, comp = c(1,2), cex = 3)

auc.splsda.MST = auroc(final.splsda.MST, roc.comp = 1, print = FALSE) 
final.splsda.MST.selectVar.1 <- selectVar(final.splsda.MST, comp = 1)$name
df.micro.feces.select <- df.micro.feces[,final.splsda.MST.selectVar.1]

# combine plots by samples
plotInd.pca.bySample.combine <- ggarrange(plotInd.pca.SBP$graph, plotInd.pca.MSA$graph, plotInd.pca.MST$graph,
                                          ncol = 3, nrow = 1)
plotInd.pca.bySample.combine

plotInd.splsda.bySample.combine <- ggarrange(plotInd.splsda.SBP$graph, plotInd.splsda.MSA$graph, plotInd.splsda.MST$graph,
                                             ncol = 3, nrow = 1)
plotInd.splsda.bySample.combine

plotInd.splsda.bySample.background.combine <- ggarrange(plotInd.splsda.SBP.background$graph, 
                                                        plotInd.splsda.MSA.background$graph, 
                                                        plotInd.splsda.MST.background$graph,
                                                        ncol = 3, nrow = 1)
plotInd.splsda.bySample.background.combine

plotInd.final.splsda.bySample.combine <- ggarrange(plotInd.final.splsda.SBP$graph, 
                                                   plotInd.final.splsda.MSA$graph, 
                                                   plotInd.final.splsda.MST$graph,
                                                   ncol = 3, nrow = 1)
plotInd.final.splsda.bySample.combine


####genus####
##### subgingival plaque ####
pca.SBP.g = pca(df.micro.plaque.g, ncomp = 10, center = TRUE, scale = T) 
plot(pca.SBP.g)
plotInd.pca.SBP.g <- plotIndiv(pca.SBP.g, group = Y.sub, ind.names = FALSE, # plot the samples projected
                             legend = FALSE, title = 'Subgingival plaque')


splsda.SBP.g <- splsda(df.micro.plaque.g, Y.sub, ncomp = 10)

plotInd.splsda.SBP.g <- plotIndiv(splsda.SBP.g, comp = 1:2, 
                                group = Y.sub, ind.names = FALSE,  # colour points by class
                                ellipse = TRUE, # include 95% confidence ellipse for each class
                                legend = FALSE, title = 'Subgingival plaque')

background.SBP.g = background.predict(splsda.SBP.g, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotInd.splsda.SBP.g.background <- plotIndiv(splsda.SBP.g, comp = 1:2,
                                           group = Y.sub, ind.names = FALSE, # colour points by class
                                           background = background.SBP.g, # include prediction background for each class
                                           legend = FALSE, title = "Subgingival plaque")

perf.splsda.SBP.g <- perf(splsda.SBP.g, validation = "Mfold", 
                        folds = 10, nrepeat = 50, cpus = 14,# use repeated cross-validation
                        progressBar = T, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.SBP.g, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.SBP.g$choice.ncomp


list.keepX.splsda.SBP.g <- c(seq(10, 45, 5), seq(50,140,10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.SBP.g <- tune.splsda(df.micro.plaque.g, Y.sub, ncomp = 2, 
                               validation = 'Mfold',
                               folds = 10, nrepeat = 50, # use repeated cross-validation
                               dist = 'max.dist', # use max.dist measure
                               measure = "BER", # use balanced error rate of dist measure
                               test.keepX = list.keepX.splsda.SBP.g,
                               progressBar = T,
                               cpus = 14) 

plot(tune.splsda.SBP.g, col = color.jet(2))

tune.splsda.SBP.g$choice.ncomp$ncomp #1
tune.splsda.SBP.g$choice.keepX

optimal.ncomp.splsda.SBP.g <- tune.splsda.SBP.g$choice.ncomp$ncomp
optimal.ncomp.splsda.SBP.g <- 2
optimal.keepX.splsda.SBP.g <- tune.splsda.SBP.g$choice.keepX[1:optimal.ncomp.splsda.SBP.g]

final.splsda.SBP.g <- splsda(df.micro.plaque.g, Y.sub, 
                           ncomp = optimal.ncomp.splsda.SBP.g, 
                           keepX = optimal.keepX.splsda.SBP.g)
biplot.final.splsda.SBP.g <- biplot(final.splsda.SBP.g, group = Y.sub, 
                                    legend.title = '') +
  theme(aspect.ratio = 1)

plotLoadings.final.splsda.SBP.g <- plotLoadings(final.splsda.SBP.g, comp = 1, contrib = 'max', method = 'median', title = "SBP")
df.pld.final.splsda.SBP.g <- plotLoadings.final.splsda.SBP.g %>%
  tibble::rownames_to_column("Genus") %>% 
  filter(GroupContrib != "tie")%>%
  left_join(taxa.micro.g, by = "Genus")
saveRDS(df.pld.final.splsda.SBP.g, "microbiome sPLSDA/df.pld.final.splsda.SBP.g.rds")
write.csv(df.pld.final.splsda.SBP.g, "microbiome sPLSDA/df.pld.final.splsda.SBP.g.csv")
pld.final.splsda.SBP.g <- ggplot(data=df.pld.final.splsda.SBP.g, 
                                 aes(x=importance, y = reorder(Genus, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Subgingival plaque\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A"), labels = c("Non-Perio", "Perio")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        title=element_text(size=16),
        axis.title = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_text(size=14, face = "italic"),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.SBP.g <- plotIndiv(final.splsda.SBP.g, comp = c(1,2), 
                                      group = Y.sub, ind.names = FALSE, 
                                      col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                      ellipse = TRUE, legend = FALSE, 
                                      title = 'Subgingival plaque')
pl.plotInd.final.splsda.SBP.g <- plotInd.final.splsda.SBP.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)



# generate the CIM, using the legend and colouring rows by each sample's class
cim.SBP.g <-cim(final.splsda.SBP.g, row.sideColors = Y.sub.color, title = "Subgingival plaque",
                comp = 1, 
                transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.SBP.g.final <- perf(final.splsda.SBP.g, 
                              folds = 10, nrepeat = 50, cpus = 14,
                              validation = "Mfold", dist = "max.dist", 
                              progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.SBP.g.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)


par(mfrow=c(1,1))

plotVar(final.splsda.SBP.g, comp = c(1,2), cex = 3)

auc.splsda.SBP.g = auroc(final.splsda.SBP.g, roc.comp = 1, print = FALSE) 
final.splsda.SBP.g.selectVar.1 <- selectVar(final.splsda.SBP.g, comp = 1)$name
df.micro.plaque.g.select <- df.micro.plaque.g[,final.splsda.SBP.g.selectVar.1]

##### saliva ####
pca.MSA.g = pca(df.micro.saliva.g, ncomp = 10, center = T, scale = T) 
plot(pca.MSA.g)
plotInd.pca.MSA.g <- plotIndiv(pca.MSA.g, group = Y.sub, ind.names = FALSE, 
                             legend = FALSE, title = 'Saliva')


splsda.MSA.g <- splsda(df.micro.saliva.g, Y.sub, ncomp = 10)

plotInd.splsda.MSA.g <- plotIndiv(splsda.MSA.g, comp = 1:2, 
                                group = Y.sub, ind.names = FALSE,  
                                ellipse = TRUE, 
                                legend = FALSE, title = 'Saliva')

background.MSA.g = background.predict(splsda.MSA.g, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MSA.g.background <- plotIndiv(splsda.MSA.g, comp = 1:2,
                                           group = Y.sub, ind.names = FALSE, 
                                           background = background.MSA.g, 
                                           legend = FALSE, title = "Saliva")

perf.splsda.MSA.g<- perf(splsda.MSA.g, validation = "Mfold", 
                       folds = 10, nrepeat = 50, cpus = 14,
                       progressBar = T, auc = TRUE) 

plot(perf.splsda.MSA.g, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MSA.g$choice.ncomp


list.keepX.splsda.MSA.g <- c(seq(10, 45, 5), seq(50,120,10))


tune.splsda.MSA.g <- tune.splsda(df.micro.saliva.g, Y.sub, ncomp = 3, 
                               validation = 'Mfold',
                               folds = 10, nrepeat = 50, 
                               dist = 'max.dist', 
                               measure = "BER", 
                               test.keepX = list.keepX.splsda.MSA.g,
                               progressBar = T,
                               cpus = 14) 

plot(tune.splsda.MSA.g, col = color.jet(3))

tune.splsda.MSA.g$choice.ncomp$ncomp #2
tune.splsda.MSA.g$choice.keepX

optimal.ncomp.splsda.MSA.g <- tune.splsda.MSA.g$choice.ncomp$ncomp
optimal.ncomp.splsda.MSA.g <- 2
optimal.keepX.splsda.MSA.g <- tune.splsda.MSA.g$choice.keepX[1:optimal.ncomp.splsda.MSA.g]

final.splsda.MSA.g <- splsda(df.micro.saliva.g, Y.sub, 
                           ncomp = optimal.ncomp.splsda.MSA.g, 
                           keepX = optimal.keepX.splsda.MSA.g)
biplot.final.splsda.MSA.g <- biplot(final.splsda.MSA.g, group = Y.sub, 
                                    legend.title = '') +
  theme(aspect.ratio = 1)

plotLoadings.final.splsda.MSA.g <- plotLoadings(final.splsda.MSA.g, comp = 1, contrib = 'max', method = 'median', title = "MSA")
df.pld.final.splsda.MSA.g <- plotLoadings.final.splsda.MSA.g %>%
  tibble::rownames_to_column("Genus") %>% 
  filter(GroupContrib != "tie")%>%
  left_join(taxa.micro.g, by = "Genus")
saveRDS(df.pld.final.splsda.MSA.g, "microbiome sPLSDA/df.pld.final.splsda.MSA.g.rds")
write.csv(df.pld.final.splsda.MSA.g, "microbiome sPLSDA/df.pld.final.splsda.MSA.g.csv")
pld.final.splsda.MSA.g <- ggplot(data=df.pld.final.splsda.MSA.g, 
                                 aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text=element_text(size=12, face = "italic"),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotLoadings.final.splsda.MSA.g.comp2 <- plotLoadings(final.splsda.MSA.g, comp = 2, contrib = 'max', method = 'median', title = "MSA")
df.pld.final.splsda.MSA.g.comp2 <- plotLoadings.final.splsda.MSA.g.comp2 %>%
  tibble::rownames_to_column("ASV") %>% 
  filter(GroupContrib != "tie")
saveRDS(df.pld.final.splsda.MSA.g.comp2, "microbiome sPLSDA/df.pld.final.splsda.MSA.g.comp2.rds")
pld.final.splsda.MSA.g.comp2 <- ggplot(data=df.pld.final.splsda.MSA.g.comp2, 
                                 aes(x=importance, y = reorder(Genus, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Saliva\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.MSA.g <- plotIndiv(final.splsda.MSA.g, comp = c(1,2), 
                                        group = Y.sub, ind.names = FALSE, 
                                        col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                        ellipse = TRUE, legend = FALSE, 
                                        title = 'Saliva')
pl.plotInd.final.splsda.MSA.g <- plotInd.final.splsda.MSA.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)



cim.MSA.g <- cim(final.splsda.MSA.g, row.sideColors = Y.sub.color, title = "Saliva", 
                 comp = 1,
                 transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.MSA.g.final <- perf(final.splsda.MSA.g, 
                              folds = 10, nrepeat = 50, cpus = 14,
                              validation = "Mfold", dist = "max.dist", 
                              progressBar = TRUE)

par(mfrow=c(1,2))
plot(perf.splsda.MSA.g.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.MSA.g.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.MSA.g, comp = c(1,2), cex = 3)

auc.splsda.MSA.g = auroc(final.splsda.MSA.g, roc.comp = 1, print = FALSE) 
final.splsda.MSA.g.selectVar.1 <- selectVar(final.splsda.MSA.g, comp = 1)$name
final.splsda.MSA.g.selectVar.2 <- selectVar(final.splsda.MSA.g, comp = 2)$name
df.micro.saliva.g.select <- df.micro.saliva.g[,unique(c(final.splsda.MSA.g.selectVar.1,
                                                 final.splsda.MSA.g.selectVar.2))]
##### feces ####

pca.MST.g = pca(df.micro.feces.g, ncomp = 10, center = T, scale = T) 
plot(pca.MST.g)
plotInd.pca.MST.g <- plotIndiv(pca.MST.g, group = Y.sub, ind.names = FALSE, 
                             legend = FALSE, title = 'Feces')


splsda.MST.g <- splsda(df.micro.feces.g, Y.sub, ncomp = 10)
plotInd.splsda.MST.g <- plotIndiv(splsda.MST.g, comp = 1:2, 
                                group = Y.sub, ind.names = FALSE,  
                                ellipse = TRUE, 
                                legend = FALSE, title = 'Feces')

background.MST.g = background.predict(splsda.MST.g, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MST.g.background <- plotIndiv(splsda.MST.g, comp = 1:2,
                                           group = Y.sub, ind.names = FALSE, 
                                           background = background.MST.g, 
                                           legend = FALSE, title = "Feces")

perf.splsda.MST.g <- perf(splsda.MST.g, validation = "Mfold", 
                        folds = 10, nrepeat = 50, cpus = 14,
                        progressBar = T, auc = TRUE) 

plot(perf.splsda.MST.g, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MST.g$choice.ncomp


list.keepX.splsda.MST.g <- c(seq(10, 45, 5), seq(50,150,10))


tune.splsda.MST.g <- tune.splsda(df.micro.feces.g, Y.sub, ncomp = 2, 
                               validation = 'Mfold',
                               folds = 10, nrepeat = 50, 
                               dist = 'max.dist', 
                               measure = "BER", 
                               test.keepX = list.keepX.splsda.MST.g,
                               progressBar = T,
                               cpus = 14) 

plot(tune.splsda.MST.g, col = color.jet(2))

tune.splsda.MST.g$choice.ncomp$ncomp #2
tune.splsda.MST.g$choice.keepX

optimal.ncomp.splsda.MST.g <- tune.splsda.MST.g$choice.ncomp$ncomp
#optimal.ncomp.splsda.MST.g <- 2
optimal.keepX.splsda.MST.g <- tune.splsda.MST.g$choice.keepX[1:optimal.ncomp.splsda.MST.g]

final.splsda.MST.g <- splsda(df.micro.feces.g, Y.sub, 
                           ncomp = optimal.ncomp.splsda.MST.g, 
                           keepX = optimal.keepX.splsda.MST.g)
biplot.final.splsda.MST.g <- biplot(final.splsda.MST.g, group = Y.sub, col.per.group = c("#00468B", "#AD002A"),
                                    legend.title = '') +
  theme(aspect.ratio = 1)

plotLoadings.final.splsda.MST.g <- plotLoadings(final.splsda.MST.g, comp = 1, contrib = 'max', method = 'median', title = "MST")
df.pld.final.splsda.MST.g <- plotLoadings.final.splsda.MST.g %>%
  tibble::rownames_to_column("Genus") %>% 
  filter(GroupContrib != "tie")%>%
  left_join(taxa.micro.g, by = "Genus")
saveRDS(df.pld.final.splsda.MST.g, "microbiome sPLSDA/df.pld.final.splsda.MST.g.rds")
write.csv(df.pld.final.splsda.MST.g, "microbiome sPLSDA/df.pld.final.splsda.MST.g.csv")
pld.final.splsda.MST.g <- ggplot(data=df.pld.final.splsda.MST.g, 
                                 aes(x=importance, y = reorder(Genus, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.3), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        legend.position = "none",
        title=element_text(size=16),
        axis.title = element_blank(),
        axis.text.x =element_text(size=14),
        axis.text.y =element_text(size=14, face = "italic"),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.MST.g <- plotIndiv(final.splsda.MST.g, comp = c(1,2), 
                                        group = Y.sub, ind.names = FALSE, 
                                        col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                        ellipse = TRUE, legend = FALSE, 
                                        title = 'Feces')
pl.plotInd.final.splsda.MST.g <- plotInd.final.splsda.MST.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)



cim.MST.g <- cim(final.splsda.MST.g, row.sideColors = Y.sub.color, title = "Feces",
                 comp = 1,
                 transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)
plotVar(final.splsda.MST.g, comp = c(1,2), cex = 3)

auc.splsda.MST.g = auroc(final.splsda.MST.g, roc.comp = 1, print = FALSE) 
final.splsda.MST.g.selectVar.1 <- selectVar(final.splsda.MST.g, comp = 1)$name
df.micro.feces.g.select <- df.micro.feces.g[,final.splsda.MST.g.selectVar.1]

####species####
##### subgingival plaque ####
pca.SBP.s = pca(df.micro.plaque.s, ncomp = 10, center = TRUE, scale = T) 
plot(pca.SBP.s)
plotInd.pca.SBP.s <- plotIndiv(pca.SBP.s, group = Y.sub, ind.names = FALSE, # plot the samples projected
                               legend = FALSE, title = 'Subgingival plaque')


splsda.SBP.s <- splsda(df.micro.plaque.s, Y.sub, ncomp = 10)

plotInd.splsda.SBP.s <- plotIndiv(splsda.SBP.s, comp = 1:2, 
                                  group = Y.sub, ind.names = FALSE,  # colour points by class
                                  ellipse = TRUE, # include 95% confidence ellipse for each class
                                  legend = FALSE, title = 'Subgingival plaque')

background.SBP.s = background.predict(splsda.SBP.s, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotInd.splsda.SBP.s.background <- plotIndiv(splsda.SBP.s, comp = 1:2,
                                             group = Y.sub, ind.names = FALSE, # colour points by class
                                             background = background.SBP.s, # include prediction background for each class
                                             legend = FALSE, title = "Subgingival plaque")

perf.splsda.SBP.s <- perf(splsda.SBP.s, validation = "Mfold", 
                          folds = 10, nrepeat = 50, cpus = 14,# use repeated cross-validation
                          progressBar = T, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.SBP.s, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.SBP.s$choice.ncomp


list.keepX.splsda.SBP.s <- c(seq(10, 45, 5), seq(50,140,10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.SBP.s <- tune.splsda(df.micro.plaque.s, Y.sub, ncomp = 2, 
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX.splsda.SBP.s,
                                 progressBar = T,
                                 cpus = 14) 

plot(tune.splsda.SBP.s, col = color.jet(2))

tune.splsda.SBP.s$choice.ncomp$ncomp #1
tune.splsda.SBP.s$choice.keepX

optimal.ncomp.splsda.SBP.s <- tune.splsda.SBP.s$choice.ncomp$ncomp
optimal.ncomp.splsda.SBP.s <- 2
optimal.keepX.splsda.SBP.s <- tune.splsda.SBP.s$choice.keepX[1:optimal.ncomp.splsda.SBP.s]

final.splsda.SBP.s <- splsda(df.micro.plaque.s, Y.sub, 
                             ncomp = optimal.ncomp.splsda.SBP.s, 
                             keepX = optimal.keepX.splsda.SBP.s)
plotLoadings(final.splsda.SBP.s, comp = 1, contrib = 'max', method = 'median', title = "SBP")

plotLoadings.final.splsda.SBP.s <- plotLoadings(final.splsda.SBP.s, comp = 1, contrib = 'max', method = 'median', title = "SBP")
df.pld.final.splsda.SBP.s <- plotLoadings.final.splsda.SBP.s %>%
  tibble::rownames_to_column("Species") %>% 
  filter(GroupContrib != "tie") %>%
  left_join(taxa.micro.sp, by = "Species")
saveRDS(df.pld.final.splsda.SBP.s, "microbiome sPLSDA/df.pld.final.splsda.SBP.s.rds")
write.csv(df.pld.final.splsda.SBP.s, "microbiome sPLSDA/df.pld.final.splsda.SBP.s.csv")
pld.final.splsda.SBP.s <- ggplot(data=df.pld.final.splsda.SBP.s, 
                                 aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Subgingival plaque\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.6, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.SBP.s <- plotIndiv(final.splsda.SBP.s, comp = c(1,2), 
                                        group = Y.sub, ind.names = FALSE, 
                                        col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                        ellipse = TRUE, legend = FALSE, 
                                        title = 'Subgingival plaque')
pl.plotInd.final.splsda.SBP.s <- plotInd.final.splsda.SBP.s$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

# generate the CIM, using the legend and colouring rows by each sample's class
cim.SBP.s <-cim(final.splsda.SBP.s, row.sideColors = Y.sub.color, title = "Subgingival plaque",
                comp = 1,
                transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.SBP.s.final <- perf(final.splsda.SBP.s, 
                                folds = 10, nrepeat = 50, cpus = 14,
                                validation = "Mfold", dist = "max.dist", 
                                progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.SBP.s.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.SBP.s.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.SBP.s, comp = c(1,2), cex = 3)

auc.splsda.SBP.s = auroc(final.splsda.SBP.s, roc.comp = 1, print = FALSE) 
final.splsda.SBP.s.selectVar.1 <- selectVar(final.splsda.SBP.s, comp = 1)$name
df.micro.plaque.s.select <- df.micro.plaque.s[,final.splsda.SBP.s.selectVar.1]

##### saliva ####
pca.MSA.s = pca(df.micro.saliva.s, ncomp = 10, center = T, scale = T) 
plot(pca.MSA.s)
plotInd.pca.MSA.s <- plotIndiv(pca.MSA.s, group = Y.sub, ind.names = FALSE, 
                               legend = FALSE, title = 'Saliva')


splsda.MSA.s <- splsda(df.micro.saliva.s, Y.sub, ncomp = 10)

plotInd.splsda.MSA.s <- plotIndiv(splsda.MSA.s, comp = 1:2, 
                                  group = Y.sub, ind.names = FALSE,  
                                  ellipse = TRUE, 
                                  legend = FALSE, title = 'Saliva')

background.MSA.s = background.predict(splsda.MSA.s, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MSA.s.background <- plotIndiv(splsda.MSA.s, comp = 1:2,
                                             group = Y.sub, ind.names = FALSE, 
                                             background = background.MSA.s, 
                                             legend = FALSE, title = "Saliva")

perf.splsda.MSA.s<- perf(splsda.MSA.s, validation = "Mfold", 
                         folds = 10, nrepeat = 50, cpus = 14,
                         progressBar = T, auc = TRUE) 

plot(perf.splsda.MSA.s, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MSA.s$choice.ncomp


list.keepX.splsda.MSA.s <- c(seq(10, 45, 5), seq(50,130,10))


tune.splsda.MSA.s <- tune.splsda(df.micro.saliva.s, Y.sub, ncomp = 3, 
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, 
                                 dist = 'max.dist', 
                                 measure = "BER", 
                                 test.keepX = list.keepX.splsda.MSA.s,
                                 progressBar = T,
                                 cpus = 14) 

plot(tune.splsda.MSA.s, col = color.jet(3))

tune.splsda.MSA.s$choice.ncomp$ncomp #2
tune.splsda.MSA.s$choice.keepX

optimal.ncomp.splsda.MSA.s <- tune.splsda.MSA.s$choice.ncomp$ncomp
optimal.ncomp.splsda.MSA.s <- 2
optimal.keepX.splsda.MSA.s <- tune.splsda.MSA.s$choice.keepX[1:optimal.ncomp.splsda.MSA.s]

final.splsda.MSA.s <- splsda(df.micro.saliva.s, Y.sub, 
                             ncomp = optimal.ncomp.splsda.MSA.s, 
                             keepX = optimal.keepX.splsda.MSA.s)
plotLoadings.final.splsda.MSA.s <- plotLoadings(final.splsda.MSA.s, comp = 1, contrib = 'max', method = 'median', title = "MSA")
df.pld.final.splsda.MSA.s <- plotLoadings.final.splsda.MSA.s %>%
  tibble::rownames_to_column("Species") %>% 
  filter(GroupContrib != "tie") %>%
  left_join(taxa.micro.sp, by = "Species")
saveRDS(df.pld.final.splsda.MSA.s, "microbiome sPLSDA/df.pld.final.splsda.MSA.s.rds")
write.csv(df.pld.final.splsda.MSA.s, "microbiome sPLSDA/df.pld.final.splsda.MSA.s.csv")
pld.final.splsda.MSA.s <- ggplot(data=df.pld.final.splsda.MSA.s, 
                                 aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
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

plotInd.final.splsda.MSA.s <- plotIndiv(final.splsda.MSA.s, comp = c(1,2), 
                                        group = Y.sub, ind.names = FALSE, 
                                        col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                        ellipse = TRUE, legend = FALSE, 
                                        title = 'Saliva')
pl.plotInd.final.splsda.MSA.s <- plotInd.final.splsda.MSA.s$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)



cim.MSA.s <- cim(final.splsda.MSA.s, row.sideColors = Y.sub.color, title = "Saliva", 
                 comp = 1,
                 transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.MSA.s.final <- perf(final.splsda.MSA.s, 
                                folds = 10, nrepeat = 50, cpus = 14,
                                validation = "Mfold", dist = "max.dist", 
                                progressBar = TRUE)

par(mfrow=c(1,2))
plot(perf.splsda.MSA.s.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.MSA.s.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.MSA.s, comp = c(1,2), cex = 3)

auc.splsda.MSA.s = auroc(final.splsda.MSA.s, roc.comp = 1, print = FALSE) 
final.splsda.MSA.s.selectVar.1 <- selectVar(final.splsda.MSA.s, comp = 1)$name
final.splsda.MSA.s.selectVar.2 <- selectVar(final.splsda.MSA.s, comp = 2)$name
df.micro.saliva.s.select <- df.micro.saliva.s[,unique(c(final.splsda.MSA.s.selectVar.1,
                                                        final.splsda.MSA.s.selectVar.2))]
##### feces ####

pca.MST.s = pca(df.micro.feces.s, ncomp = 10, center = T, scale = T) 
plot(pca.MST.s)
plotInd.pca.MST.s <- plotIndiv(pca.MST.s, group = Y.sub, ind.names = FALSE, 
                               legend = FALSE, title = 'Feces')


splsda.MST.s <- splsda(df.micro.feces.s, Y.sub, ncomp = 10)
plotInd.splsda.MST.s <- plotIndiv(splsda.MST.s, comp = 1:2, 
                                  group = Y.sub, ind.names = FALSE,  
                                  ellipse = TRUE, 
                                  legend = FALSE, title = 'Feces')

background.MST.s = background.predict(splsda.MST.s, comp.predicted=2, dist = "max.dist")


plotInd.splsda.MST.s.background <- plotIndiv(splsda.MST.s, comp = 1:2,
                                             group = Y.sub, ind.names = FALSE, 
                                             background = background.MST.s, 
                                             legend = FALSE, title = "Feces")

perf.splsda.MST.s <- perf(splsda.MST.s, validation = "Mfold", 
                          folds = 10, nrepeat = 50, cpus = 14,
                          progressBar = T, auc = TRUE) 

plot(perf.splsda.MST.s, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.MST.s$choice.ncomp


list.keepX.splsda.MST.s <- c(seq(10, 45, 5), seq(50,120,10))


tune.splsda.MST.s <- tune.splsda(df.micro.feces.s, Y.sub, ncomp = 2, 
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, 
                                 dist = 'max.dist', 
                                 measure = "BER", 
                                 test.keepX = list.keepX.splsda.MST.s,
                                 progressBar = T,
                                 cpus = 14) 

plot(tune.splsda.MST.s, col = color.jet(2))

tune.splsda.MST.s$choice.ncomp$ncomp #2
tune.splsda.MST.s$choice.keepX

optimal.ncomp.splsda.MST.s <- tune.splsda.MST.s$choice.ncomp$ncomp
optimal.ncomp.splsda.MST.s <- 2
optimal.keepX.splsda.MST.s <- tune.splsda.MST.s$choice.keepX[1:optimal.ncomp.splsda.MST.s]

final.splsda.MST.s <- splsda(df.micro.feces.s, Y.sub, 
                             ncomp = optimal.ncomp.splsda.MST.s, 
                             keepX = optimal.keepX.splsda.MST.s)
plotLoadings.final.splsda.MST.s <- plotLoadings(final.splsda.MST.s, comp = 1, contrib = 'max', method = 'median', title = "MST")
df.pld.final.splsda.MST.s <- plotLoadings.final.splsda.MST.s %>%
  tibble::rownames_to_column("Species") %>% 
  filter(GroupContrib != "tie") %>%
  left_join(taxa.micro.sp, by = "Species")
saveRDS(df.pld.final.splsda.MST.s, "microbiome sPLSDA/df.pld.final.splsda.MST.s.rds")
write.csv(df.pld.final.splsda.MST.s, "microbiome sPLSDA/df.pld.final.splsda.MST.s.csv")
pld.final.splsda.MST.s <- ggplot(data=df.pld.final.splsda.MST.s, 
                                 aes(x=importance, y = reorder(ASV, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="Feces\nmicrobiome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.3), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plotInd.final.splsda.MST.s <- plotIndiv(final.splsda.MST.s, comp = c(1,2), 
                                        group = Y.sub, ind.names = FALSE, 
                                        col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                        ellipse = TRUE, legend = FALSE, 
                                        title = 'Saliva')
pl.plotInd.final.splsda.MST.s <- plotInd.final.splsda.MST.s$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)



cim.MST.s <- cim(final.splsda.MST.s, row.sideColors = Y.sub.color, title = "Feces",
                 comp = 1,
                 transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.MST.s.final <- perf(final.splsda.MST.s, 
                                folds = 10, nrepeat = 50, cpus = 14,
                                validation = "Mfold", dist = "max.dist", 
                                progressBar = TRUE)

par(mfrow=c(1,2))
plot(perf.splsda.MST.s.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.MST.s.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

plotVar(final.splsda.MST.s, comp = c(1,2), cex = 3)

auc.splsda.MST.s = auroc(final.splsda.MST.s, roc.comp = 1, print = FALSE) 
final.splsda.MST.s.selectVar.1 <- selectVar(final.splsda.MST.s, comp = 1)$name
final.splsda.MST.s.selectVar.2 <- selectVar(final.splsda.MST.s, comp = 2)$name
df.micro.feces.s.select <- df.micro.feces.s[,unique(c(final.splsda.MST.s.selectVar.1,
                                                        final.splsda.MST.s.selectVar.2))]


####AUC Cliff####
library(ROCR)
library(rcompanion)
library(ggpubr)
#####genus####
#MST
roc.perf.feces.g <- list()
auc.perf.feces.g <- list()
df_for_ROC_Cliff.g <- data.frame(df.micro.feces.g.select, Disease = Y.sub)
df.auc.cliff.feces.g <- data.frame(matrix(ncol=2,nrow=ncol(df.micro.feces.g.select), 
                                  dimnames=list(colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)]), 
                                                c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)])){
  pred <- prediction(df_for_ROC_Cliff.g[,i],df_for_ROC_Cliff.g$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.feces.g[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.feces.g[[i]] <- auc.perf
  df.auc.cliff.feces.g[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.feces.g[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                               data = df_for_ROC_Cliff.g)
}
rownames(df.auc.cliff.feces.g) <- colnames(df.micro.feces.g.select)
df.auc.cliff.feces.GroupContrib.g <- merge(df.auc.cliff.feces.g, 
                                                    plotLoadings.final.splsda.MST.g, by = 0)
df.auc.cliff.feces.GroupContrib.g$GroupContrib <- factor(df.auc.cliff.feces.GroupContrib.g$GroupContrib, 
                                                       levels = c("NPD", "PD", "tie"))
my_comparisons <- list( c("NPD", "PD"))
shapiro.test(df.auc.cliff.feces.GroupContrib.g$AUC)
shapiro.test(df.auc.cliff.feces.GroupContrib.g$Cliff)

pl.AUC.feces.g <- ggboxplot(df.auc.cliff.feces.GroupContrib.g, 
          x="GroupContrib", y="AUC", color = "GroupContrib",
          fill = "GroupContrib", alpha = 0.2,
          palette = c("#00468B", "#AD002A", "#42B540"),
          notch = T, add = "jitter") + 
  labs(title = "Feces", x ="Group of Contribution", color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.8)
pl.Cliff.feces.g <- ggboxplot(df.auc.cliff.feces.GroupContrib.g, 
                          x="GroupContrib", y="Cliff", color = "GroupContrib",
                          fill = "GroupContrib", alpha = 0.2,
                          palette = c("#00468B", "#AD002A", "#42B540"),
                          notch = T, add = "jitter") + 
  labs(title = "Feces", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.6)

#SBP
roc.perf.plaque.g <- list()
auc.perf.plaque.g <- list()
df_for_ROC_Cliff.g <- data.frame(df.micro.plaque.g.select, Disease = Y.sub)
df.auc.cliff.plaque.g <- data.frame(matrix(ncol=2,nrow=ncol(df.micro.plaque.g.select), 
                                        dimnames=list(colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)]), 
                                                      c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)])){
  pred <- prediction(df_for_ROC_Cliff.g[,i],df_for_ROC_Cliff.g$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.plaque.g[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.plaque.g[[i]] <- auc.perf
  df.auc.cliff.plaque.g[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.plaque.g[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                               data = df_for_ROC_Cliff.g)
}
rownames(df.auc.cliff.plaque.g) <- colnames(df.micro.plaque.g.select)
df.auc.cliff.plaque.GroupContrib.g <- merge(df.auc.cliff.plaque.g, 
                                         plotLoadings.final.splsda.SBP.g, by = 0)
df.auc.cliff.plaque.GroupContrib.g$GroupContrib <- factor(df.auc.cliff.plaque.GroupContrib.g$GroupContrib, 
                                                       levels = c("NPD", "PD", "tie"))

shapiro.test(df.auc.cliff.plaque.GroupContrib.g$AUC)
shapiro.test(df.auc.cliff.plaque.GroupContrib.g$Cliff)

pl.AUC.plaque.g <- ggboxplot(df.auc.cliff.plaque.GroupContrib.g, 
                          x="GroupContrib", y="AUC", color = "GroupContrib",
                          fill = "GroupContrib", alpha = 0.2,
                          palette = c("#00468B", "#AD002A", "#42B540"),
                          notch = T, add = "jitter") + 
  labs(title = "Subgingival plaque", x ="Group of Contribution", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.8)
pl.Cliff.plaque.g <- ggboxplot(df.auc.cliff.plaque.GroupContrib.g, 
                            x="GroupContrib", y="Cliff", color = "GroupContrib",
                            fill = "GroupContrib", alpha = 0.2,
                            palette = c("#00468B", "#AD002A", "#42B540"),
                            notch = T, add = "jitter") + 
  labs(title = "Subgingival plaque", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.6)

#MSA
roc.perf.saliva.g <- list()
auc.perf.saliva.g <- list()
df_for_ROC_Cliff.g<- data.frame(df.micro.saliva.g.select, Disease = Y.sub)
df.auc.cliff.saliva.g <- data.frame(matrix(ncol=2,nrow=ncol(df.micro.saliva.g.select), 
                                         dimnames=list(colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)]), 
                                                       c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff.g[1:(ncol(df_for_ROC_Cliff.g)-1)])){
  pred <- prediction(df_for_ROC_Cliff.g[,i],df_for_ROC_Cliff.g$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.saliva.g[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.saliva.g[[i]] <- auc.perf
  df.auc.cliff.saliva.g[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.saliva.g[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                                data = df_for_ROC_Cliff.g)
}
rownames(df.auc.cliff.saliva.g) <- colnames(df.micro.saliva.g.select)
df.auc.cliff.saliva.GroupContrib.g <- merge(df.auc.cliff.saliva.g, 
                                          plotLoadings.final.splsda.MSA.g, by = 0)
df.auc.cliff.saliva.GroupContrib.g$GroupContrib <- factor(df.auc.cliff.saliva.GroupContrib.g$GroupContrib, 
                                                        levels = c("NPD", "PD", "tie"))

shapiro.test(df.auc.cliff.saliva.GroupContrib.g$AUC)
shapiro.test(df.auc.cliff.saliva.GroupContrib.g$Cliff)

pl.AUC.saliva.g <- ggboxplot(df.auc.cliff.saliva.GroupContrib.g, 
                           x="GroupContrib", y="AUC", color = "GroupContrib",
                           fill = "GroupContrib", alpha = 0.2,
                           palette = c("#00468B", "#AD002A", "#42B540"),
                           notch = T, add = "jitter") + 
  labs(title = "Saliva", x ="Group of Contribution", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.8)
pl.Cliff.saliva.g <- ggboxplot(df.auc.cliff.saliva.GroupContrib.g, 
                             x="GroupContrib", y="Cliff", color = "GroupContrib",
                             fill = "GroupContrib", alpha = 0.2, 
                             palette = c("#00468B", "#AD002A", "#42B540"),
                             notch = T, add = "jitter") + 
  labs(title = "Saliva", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  ylim(-0.6,0.6) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  stat_compare_means(method = "anova", label.y = 0.6)




ggarrange(ggarrange(pld.final.splsda.SBP.g,
          pld.final.splsda.MST.g,
          ncol = 2, nrow = 1,
          common.legend = TRUE),
          ggarrange(pl.AUC.plaque.g + theme(legend.position='hidden'), 
                    pl.AUC.saliva.g + theme(legend.position='hidden'), 
                    pl.AUC.feces.g + theme(legend.position='hidden'), 
                    pl.Cliff.plaque.g + theme(legend.position='hidden'), 
                    pl.Cliff.saliva.g + theme(legend.position='hidden'), 
                    pl.Cliff.feces.g + theme(legend.position='hidden'), 
                    ncol = 3, nrow = 2, legend = "bottom",
                    common.legend = TRUE),
          ncol = 1, nrow = 2)





