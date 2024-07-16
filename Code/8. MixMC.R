library(mixOmics) # import the mixOmics library
library(dplyr)
library(readxl)
library(mia)
library(stringr)
library(ggplot2)
library(igraph)

setwd("")
set.seed(123)
tse.micro <- readRDS("")


##MixMC
X.samples <- tse.micro %>%
  subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05) %>% #filter with 0.05% detection & 5% prevalence
  transformSamples(method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>% #transform counts to rubust clr
  assay("rclr") %>%
  t()
Y.samples <- sample_meta$Sample %>%
  factor(levels = c("SBP", "MSA", "MST"))
levels(Y.samples) <- c("Subgingival plaques", "Saliva", "Feces")
sample_meta <- data.frame(colData(tse.micro))


sample.pca = pca(X.samples, ncomp = 10) # undergo PCA with 10 comps
plot(sample.pca) # plot explained variance

plotIndiv(sample.pca, 
          ind.names = FALSE, # not showing sample names
          group = Y.samples, # color according to Y
          legend = TRUE,
          title = 'Sample types, PCA Comps 1&2')

basic.sample.plsda = plsda(X.samples, Y.samples,  
                           ncomp = nlevels(Y.samples))
# assess the performance of the sPLS-DA model using repeated CV
basic.sample.perf.plsda = perf(basic.sample.plsda,  
                               validation = 'Mfold', 
                               folds = 10, nrepeat = 10, 
                               progressBar = T)

# extract the optimal component number
optimal.ncomp.sample <- basic.sample.perf.plsda$choice.ncomp["BER", "max.dist"] 
plot(basic.sample.perf.plsda, overlay = 'measure', sd=TRUE) 


grid.keepX.sample = c(seq(5,500, 25))

sample.tune.splsda = tune.splsda(X.samples, Y.samples,
                                 ncomp = optimal.ncomp.sample, # use optimal component number
                                 test.keepX = grid.keepX.sample,
                                 validation = c('Mfold'),
                                 folds = 10, nrepeat = 10, # use repeated CV
                                 dist = 'max.dist', # maximum distance as metric
                                 cpus = 10,
                                 progressBar = T)

# extract the optimal component number and optimal feature count per component
optimal.keepX.sample = sample.tune.splsda$choice.keepX[1:2] 
optimal.ncomp.sample = sample.tune.splsda$choice.ncomp$ncomp 

plot(sample.tune.splsda) # plot this tuning
#final model
sample.splsda = splsda(X.samples, Y.samples, # form final sPLS-DA model
                       ncomp = optimal.ncomp.sample, 
                       keepX = optimal.keepX.sample)

plot.sample.splsda <- plotIndiv(sample.splsda,
                                comp = c(1,2),
                                ind.names = FALSE,
                                ellipse = TRUE, # include confidence ellipses
                                abline = T,
                                col = c("#0f77c1", "#00a087", "#e64b35"),
                                legend = TRUE,
                                legend.title = "",
                                title = '')
pl.sample.splsda <- plot.sample.splsda$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

legend.sample = list(legend = levels(Y.samples), # set of classes
                     col = c("#0f77c1", "#00a087", "#e64b35"), # set of colours
                     cex = 0.7) # legend size
color.sample <- Y.samples
levels(color.sample) <- c("#0f77c1", "#00a087", "#e64b35")
cim(sample.splsda,
    comp = c(1,2), margins = c(2, 20),
    row.sideColors = color.sample, row.cex = 1, 
    legend = legend.sample, 
    row.names = F, col.names = T, transpose = T,
    title = '')

plotVar(sample.splsda,
        comp = c(1,2),
        cutoff = 0.7, rad.in = 0.7,
        title = 'Sample types, Correlation Circle Plot Comps 1&2')

sample.perf.splsda = perf(sample.splsda, validation = 'Mfold', 
                          folds = 10, nrepeat = 10, cpus = 10,
                          progressBar = T, dist = 'max.dist') 
sample.perf.splsda$error.rate



plotLoadings(sample.splsda, comp = 1, 
             method = 'mean', contrib = 'max',  
             size.name = 0.7, legend = FALSE,  
             ndisplay = 40, 
             title = "(a) Loadings of first component")

plotLoadings(sample.splsda, comp = 2, 
             method = 'mean', contrib = 'max',   
             size.name = 0.7,
             ndisplay = 40, max.name.length = 50,
             title = "(b) Loadings of second comp.")

sample.splsda.selectVar.1 <- selectVar(sample.splsda, comp = 1)$name
sample.splsda.selectVar.2 <- selectVar(sample.splsda, comp = 2)$name
sample.splsda.selectVar <- c(sample.splsda.selectVar.1, sample.splsda.selectVar.2)