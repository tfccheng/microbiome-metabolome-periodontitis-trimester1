## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library
library(tidyverse)
library(vegan)
setwd("")
df_metab_feces.X.log.pareto <- readRDS("Data prepare/df_metab_feces_X.log.pareto.rds")
df_metab_saliva.X.log.pareto <- readRDS("Data prepare/df_metab_saliva_X.log.pareto.rds")
df_metab_serum.X.log.pareto <- readRDS("Data prepare/df_metab_serum_X.log.pareto.rds")

df_metab_feces_sig.X.log.pareto <- readRDS("Data prepare/df_metab_feces_X.log.pareto.sig.rds")
df_metab_saliva_sig.X.log.pareto <- readRDS("Data prepare/df_metab_saliva_X.log.pareto.sig.rds")
df_metab_serum_sig.X.log.pareto <- readRDS("Data prepare/df_metab_serum_X.log.pareto.sig.rds")

# load metabolite identification
df_metab_feces_identification <- readRDS("Data filtering/df_metab_feces_identification.rds")
df_metab_saliva_identification <- readRDS("Data filtering/df_metab_saliva_identification.rds")
df_metab_serum_identification <- readRDS("Data filtering/df_metab_serum_identification.rds")

#load meta

dim(df_metab_feces.X.log.pareto) 
meta <- readRDS("Data filtering/meta.metab.micro.combine.saliva.rds")
Y.sub <- meta$Disease
summary(Y.sub)
Y.sub.color <- case_when(Y.sub == "NPD" ~ "#00468B", 
                         Y.sub == "PD" ~ "#AD002A")
df.clinical.num.impute <- readRDS("Data prepare/df.clinical.num.impute.rds")
####GLM####
df_metab_feces.glm <- df_metab_feces.X.log.pareto %>%
  tibble::add_column(Disease = meta$Disease, .before = 1) %>%
  tibble::add_column(BMI = df.clinical.num.impute$BMI, .after = 1)
glm_metab_feces <- data.frame(ID = NA,p= NA,coef=NA,p_BMI= NA,coef_BMI=NA)
for (i in 1:(ncol(df_metab_feces.glm)-2)) {
  tmp=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_feces.glm)[i+2], "`")))
  tmp.BMI=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_feces.glm)[i+2], "`"), "+ BMI"))
  glm_test <- summary(glm(tmp, family = binomial(link = "logit"), data = df_metab_feces.glm))
  glm_test.BMI <- summary(glm(tmp.BMI, family = binomial(link = "logit"), data = df_metab_feces.glm))
  glm_metab_feces[i,1] <- colnames(df_metab_feces.glm)[i+2] #variable name
  glm_metab_feces[i,2] <- glm_test$coefficients[2,4] #p value
  glm_metab_feces[i,3] <- glm_test$coefficients[2,1] #coefficient
  glm_metab_feces[i,4] <- glm_test.BMI$coefficients[2,4] #p value
  glm_metab_feces[i,5] <- glm_test.BMI$coefficients[2,1] #coefficient
}
glm_metab_feces.sig <- glm_metab_feces %>% dplyr::filter(p<0.05&p_BMI<0.05) %>%
  dplyr::left_join(df_metab_feces_identification, by = "ID") %>% 
  arrange(desc(IsMS2ident), desc(coef))
write.csv(glm_metab_feces.sig, "glm/glm_metab_feces.csv", row.names = F)

df_metab_serum.glm <- df_metab_serum.X.log.pareto %>%
  tibble::add_column(Disease = meta$Disease, .before = 1) %>%
  tibble::add_column(BMI = df.clinical.num.impute$BMI, .after = 1)
glm_metab_serum <- data.frame(ID = NA,p= NA,coef=NA,p_BMI= NA,coef_BMI=NA)
for (i in 1:(ncol(df_metab_serum.glm)-2)) {
  tmp=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_serum.glm)[i+2], "`")))
  tmp.BMI=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_serum.glm)[i+2], "`"), "+ BMI"))
  glm_test <- summary(glm(tmp, family = binomial(link = "logit"), data = df_metab_serum.glm))
  glm_test.BMI <- summary(glm(tmp.BMI, family = binomial(link = "logit"), data = df_metab_serum.glm))
  glm_metab_serum[i,1] <- colnames(df_metab_serum.glm)[i+2] #variable name
  glm_metab_serum[i,2] <- glm_test$coefficients[2,4] #p value
  glm_metab_serum[i,3] <- glm_test$coefficients[2,1] #coefficient
  glm_metab_serum[i,4] <- glm_test.BMI$coefficients[2,4] #p value
  glm_metab_serum[i,5] <- glm_test.BMI$coefficients[2,1] #coefficient
}
glm_metab_serum.sig <- glm_metab_serum %>% dplyr::filter(p<0.05&p_BMI<0.05)%>%
  dplyr::left_join(df_metab_serum_identification, by = "ID") %>% 
  arrange(desc(IsMS2ident), desc(coef))
write.csv(glm_metab_serum.sig, "glm/glm_metab_serum.csv", row.names = F)

df_metab_saliva.glm <- df_metab_saliva.X.log.pareto %>%
  tibble::add_column(Disease = meta$Disease, .before = 1) %>%
  tibble::add_column(BMI = df.clinical.num.impute$BMI, .after = 1)
glm_metab_saliva <- data.frame(ID = NA,p= NA,coef=NA,p_BMI= NA,coef_BMI=NA)
for (i in 1:(ncol(df_metab_saliva.glm)-2)) {
  tmp=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_saliva.glm)[i+2], "`")))
  tmp.BMI=as.formula(paste0("Disease ~", paste0("`", colnames(df_metab_saliva.glm)[i+2], "`"), "+ BMI"))
  glm_test <- summary(glm(tmp, family = binomial(link = "logit"), data = df_metab_saliva.glm))
  glm_test.BMI <- summary(glm(tmp.BMI, family = binomial(link = "logit"), data = df_metab_saliva.glm))
  glm_metab_saliva[i,1] <- colnames(df_metab_saliva.glm)[i+2] #variable name
  glm_metab_saliva[i,2] <- glm_test$coefficients[2,4] #p value
  glm_metab_saliva[i,3] <- glm_test$coefficients[2,1] #coefficient
  glm_metab_saliva[i,4] <- glm_test.BMI$coefficients[2,4] #p value
  glm_metab_saliva[i,5] <- glm_test.BMI$coefficients[2,1] #coefficient
}
glm_metab_saliva.sig <- glm_metab_saliva %>% dplyr::filter(p<0.05&p_BMI<0.05)%>%
  dplyr::left_join(df_metab_saliva_identification, by = "ID") %>% 
  arrange(desc(IsMS2ident), desc(coef))
write.csv(glm_metab_saliva.sig, "glm/glm_metab_saliva.csv", row.names = F)

set.seed(1234)
##PERMANOVA subsample by disease euclidean
permanova.feces <- adonis2(df_metab_feces.X.log.pareto ~ Disease,
                           by = "margin", 
                           data = meta,
                           method = "euclidean",
                           permutations = 999)

permanova.saliva <- adonis2(df_metab_saliva.X.log.pareto ~ Disease,
                                      by = "margin", 
                                      data = meta,
                                      method = "euclidean",
                                      permutations = 999)

permanova.serum <- adonis2(df_metab_serum.X.log.pareto ~ Disease,
                            by = "margin", 
                            data = meta,
                            method = "euclidean",
                            permutations = 999)

###################################################################
#######################     SPLSDA     ############################
###################################################################
##### feces ####
pca.metab.feces = pca(df_metab_feces.X.log.pareto, ncomp = 10, center = T, scale = T) 
plot(pca.metab.feces)
plotInd.pca.metab.feces <-plotIndiv(pca.metab.feces, group = Y.sub, ind.names = FALSE, ellipse = T, # plot the samples projected
                                    col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                    legend = T, title = '')
pl.plotInd.pca.metab.feces <- plotInd.pca.metab.feces$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
pca.metab.feces_sig = pca(df_metab_feces_sig.X.log.pareto, ncomp = 10, center = TRUE, scale = T) 
plot(pca.metab.feces_sig)
plotIndiv(pca.metab.feces_sig, group = Y.sub, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on fecal metabolome, comp 1 - 2')

splsda.metab.feces <- splsda(df_metab_feces.X.log.pareto, Y.sub, ncomp = 10)

plotIndiv(splsda.metab.feces, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

splsda.metab.feces_sig <- splsda(df_metab_feces_sig.X.log.pareto, Y.sub, ncomp = 10)

plotIndiv(splsda.metab.feces_sig, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')


background.metab.feces = background.predict(splsda.metab.feces, comp.predicted=2, dist = "max.dist")
background.metab.feces_sig = background.predict(splsda.metab.feces_sig, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(splsda.metab.feces, comp = 1:2,
          group = Y.sub, ind.names = FALSE, # colour points by class
          background = background.metab.feces, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")
plotIndiv(splsda.metab.feces_sig, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.metab.feces_sig, 
          legend = TRUE, title = " (b) PLSDA with prediction background")

perf.splsda.metab.feces <- perf(splsda.metab.feces, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,# use repeated cross-validation
                                progressBar = T, auc = TRUE) # include AUC values
perf.splsda.metab.feces_sig <- perf(splsda.metab.feces_sig, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,# use repeated cross-validation
                                progressBar = T, auc = TRUE)

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.metab.feces, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.feces$choice.ncomp
plot(perf.splsda.metab.feces_sig, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.feces_sig$choice.ncomp

list.keepX.splsda.metab.feces <- c(seq(10, 100, 10), seq(150,250,50), seq(300, 1000, 100))
list.keepX.splsda.metab.feces_sig <- c(seq(10, 100, 10), seq(150,250,50), c(300, 400))
# undergo the tuning process to determine the optimal number of variables
tune.splsda.metab.feces <- tune.splsda(df_metab_feces.X.log.pareto, Y.sub, ncomp = 3, 
                                       validation = 'Mfold',
                                       folds = 10, nrepeat = 50, # use repeated cross-validation
                                       dist = 'max.dist', # use max.dist measure
                                       measure = "BER", # use balanced error rate of dist measure
                                       test.keepX = list.keepX.splsda.metab.feces,
                                       progressBar = T,
                                       cpus = 14) 
tune.splsda.metab.feces_sig <- tune.splsda(df_metab_feces_sig.X.log.pareto, Y.sub, ncomp = 3, 
                                       validation = 'Mfold',
                                       folds = 10, nrepeat = 50, # use repeated cross-validation
                                       dist = 'max.dist', # use max.dist measure
                                       measure = "BER", # use balanced error rate of dist measure
                                       test.keepX = list.keepX.splsda.metab.feces_sig,
                                       progressBar = T,
                                       cpus = 14) 

plot(tune.splsda.metab.feces, col = color.jet(3))
plot(tune.splsda.metab.feces_sig, col = color.jet(3))

tune.splsda.metab.feces$choice.ncomp$ncomp #1
tune.splsda.metab.feces$choice.keepX

optimal.ncomp.splsda.metab.feces <- tune.splsda.metab.feces$choice.ncomp$ncomp #1
optimal.ncomp.splsda.metab.feces <- 2
optimal.keepX.splsda.metab.feces <- tune.splsda.metab.feces$choice.keepX[1:optimal.ncomp.splsda.metab.feces]

optimal.ncomp.splsda.metab.feces_sig <- tune.splsda.metab.feces_sig$choice.ncomp$ncomp
optimal.keepX.splsda.metab.feces_sig <- tune.splsda.metab.feces_sig$choice.keepX[1:optimal.ncomp.splsda.metab.feces_sig]

final.splsda.metab.feces <- splsda(df_metab_feces.X.log.pareto, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.metab.feces, 
                                   keepX = optimal.keepX.splsda.metab.feces)

final.splsda.metab.feces.sig <- splsda(df_metab_feces_sig.X.log.pareto, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.metab.feces_sig, 
                                   keepX = optimal.keepX.splsda.metab.feces_sig)


plotInd.final.splsda.metab.feces <- plotIndiv(final.splsda.metab.feces, comp = c(1,2), 
                                      group = Y.sub, ind.names = FALSE, 
                                      col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                      ellipse = TRUE, legend = FALSE, 
                                      title = 'Feces')
pl.plotInd.final.splsda.metab.feces <- plotInd.final.splsda.metab.feces$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotIndiv(final.splsda.metab.feces.sig, comp = c(1,2), 
          group = Y.sub, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, 
          title = ' (a) sPLS-DA on fecal metabolome, comp 1 & 2')


legend=list(legend = levels(Y.sub), # set of classes
            col = c("#00468B", "#AD002A"), # set of colours
            title = "", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim(final.splsda.metab.feces, row.sideColors = Y.sub.color, comp = 1,
    transpose = TRUE, col.names = TRUE, row.names = FALSE, legend = legend)
cim(final.splsda.metab.feces.sig, row.sideColors = Y.sub.color, comp = 1,
    transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.metab.feces.final <- perf(final.splsda.metab.feces, 
                          folds = 10, nrepeat = 50, cpus = 14,
                          validation = "Mfold", dist = "max.dist", 
                          progressBar = TRUE)

perf.splsda.metab.feces.final.sig <- perf(final.splsda.metab.feces.sig, 
                                      folds = 10, nrepeat = 50, cpus = 14,
                                      validation = "Mfold", dist = "max.dist", 
                                      progressBar = TRUE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,2))
plot(perf.splsda.metab.feces.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.feces.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(perf.splsda.metab.feces.final.sig$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.feces.final.sig$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
par(mfrow=c(1,1))

plotVar(final.splsda.metab.feces, comp = c(1,2), cex = 3)
plotVar(final.splsda.metab.feces.sig, comp = c(1,2), cex = 3)

auc.splsda.metab.feces = auroc(final.splsda.metab.feces, roc.comp = 1, print = FALSE) 
auc.splsda.metab.feces.sig = auroc(final.splsda.metab.feces.sig, roc.comp = 1, print = FALSE) 

splsda.metab.feces.selectVar.1 <- selectVar(final.splsda.metab.feces, comp = 1)$name
df_metab_feces.X.log.pareto.select <- df_metab_feces.X.log.pareto[,splsda.metab.feces.selectVar.1]


#comp1
plotLoadings.final.splsda.metab.feces <- plotLoadings(final.splsda.metab.feces, comp = 1, contrib = 'max', method = 'median', title = "feces")
df.pld.splsda.metab.feces <- plotLoadings.final.splsda.metab.feces %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_feces_identification, by = "ID") 
write.csv(df.pld.splsda.metab.feces, "metabolome sPLSDA/df.pld.splsda.metab.feces.csv", row.names = F)
df.pld.splsda.metab.feces$identification_pld <- ifelse(is.na(df.pld.splsda.metab.feces$identification), 
                                                         df.pld.splsda.metab.feces$ID, 
                                                         df.pld.splsda.metab.feces$identification)

df.pld.splsda.metab.feces.identified <- df.pld.splsda.metab.feces %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.feces.identified$identification_id <- str_trunc(df.pld.splsda.metab.feces.identified$identification_id, 100)

pld.splsda.metab.feces.identified <- ggplot(data=df.pld.splsda.metab.feces.identified, 
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


df.pld.splsda.metab.feces.MS2 <- df.pld.splsda.metab.feces %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.feces.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                           df.pld.splsda.metab.feces.MS2$identification_MS2)
write.csv(df.pld.splsda.metab.feces.MS2, "metabolome sPLSDA/df.pld.splsda.metab.feces.MS2.csv", row.names = F)
pld.splsda.metab.feces.MS2 <- ggplot(data=df.pld.splsda.metab.feces.MS2,# %>% filter(abs(importance)>0.05), 
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

df.pld.splsda.metab.feces.single.id <- df.pld.splsda.metab.feces %>%
  filter(!is.na(identification)) %>%#filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.feces.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                       df.pld.splsda.metab.feces.single.id$identification_single_id)
pld.splsda.metab.feces.single.id <- ggplot(data=df.pld.splsda.metab.feces.single.id %>% filter(abs(importance) > 0.025), 
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



##### saliva ####
pca.metab.saliva = pca(df_metab_saliva.X.log.pareto, ncomp = 10, center = TRUE, scale = T) 
plot(pca.metab.saliva)
plotIndiv(pca.metab.saliva, group = Y.sub, ind.names = FALSE, ellipse = T,
          legend = TRUE, title = 'PCA on saliva metabolome, comp 1 - 2')

pca.metab.saliva.sig = pca(df_metab_saliva_sig.X.log.pareto, ncomp = 10, center = TRUE, scale = T) 
plot(pca.metab.saliva.sig)
plotIndiv(pca.metab.saliva.sig, group = Y.sub, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on saliva metabolome, comp 1 - 2')

splsda.metab.saliva <- splsda(df_metab_saliva.X.log.pareto, Y.sub, ncomp = 10)
splsda.metab.saliva.sig <- splsda(df_metab_saliva_sig.X.log.pareto, Y.sub, ncomp = 10)

plotIndiv(splsda.metab.saliva, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

plotIndiv(splsda.metab.saliva.sig, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

background.metab.saliva = background.predict(splsda.metab.saliva, comp.predicted=2, dist = "max.dist")
background.metab.saliva.sig = background.predict(splsda.metab.saliva.sig, comp.predicted=2, dist = "max.dist")

plotIndiv(splsda.metab.saliva, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.metab.saliva, 
          legend = TRUE, title = " (b) PLSDA with prediction background")

plotIndiv(splsda.metab.saliva.sig, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.metab.saliva.sig, 
          legend = TRUE, title = " (b) PLSDA with prediction background")

perf.splsda.metab.saliva<- perf(splsda.metab.saliva, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

perf.splsda.metab.saliva.sig <- perf(splsda.metab.saliva.sig, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

plot(perf.splsda.metab.saliva, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.saliva$choice.ncomp

plot(perf.splsda.metab.saliva.sig, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.saliva.sig$choice.ncomp

list.keepX.splsda.metab.saliva <- c(seq(10, 100, 10), seq(150,250,50), seq(300, 500, 100))
list.keepX.splsda.metab.saliva.sig <- c(seq(10, 100, 10), seq(150,250,50))

tune.splsda.metab.saliva <- tune.splsda(df_metab_saliva.X.log.pareto, Y.sub, ncomp = 3, 
                                        validation = 'Mfold',
                                        folds = 20, nrepeat = 50, 
                                        dist = 'max.dist', 
                                        measure = "BER", 
                                        test.keepX = list.keepX.splsda.metab.saliva,
                                        progressBar = T,
                                        cpus = 14) 

plot(tune.splsda.metab.saliva, col = color.jet(3))

tune.splsda.metab.saliva.sig <- tune.splsda(df_metab_saliva_sig.X.log.pareto, Y.sub, ncomp = 3, 
                                        validation = 'Mfold',
                                        folds = 20, nrepeat = 50, 
                                        dist = 'max.dist', 
                                        measure = "BER", 
                                        test.keepX = list.keepX.splsda.metab.saliva.sig,
                                        progressBar = T,
                                        cpus = 14) 

plot(tune.splsda.metab.saliva.sig, col = color.jet(3))



tune.splsda.metab.saliva$choice.ncomp$ncomp
tune.splsda.metab.saliva$choice.keepX

optimal.ncomp.splsda.metab.saliva <- tune.splsda.metab.saliva$choice.ncomp$ncomp #3
#optimal.ncomp.splsda.metab.saliva <- 2
optimal.keepX.splsda.metab.saliva <- tune.splsda.metab.saliva$choice.keepX[1:optimal.ncomp.splsda.metab.saliva]

optimal.ncomp.splsda.metab.saliva.sig <- tune.splsda.metab.saliva.sig$choice.ncomp$ncomp #3
optimal.keepX.splsda.metab.saliva.sig <- tune.splsda.metab.saliva.sig$choice.keepX[1:optimal.ncomp.splsda.metab.saliva.sig]


final.splsda.metab.saliva <- splsda(df_metab_saliva.X.log.pareto, Y.sub, 
                                    ncomp = optimal.ncomp.splsda.metab.saliva, 
                                    keepX = optimal.keepX.splsda.metab.saliva)

plotInd.final.splsda.metab.saliva <- plotIndiv(final.splsda.metab.saliva, comp = c(1,2), 
                                              group = Y.sub, ind.names = FALSE, 
                                              col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                              ellipse = TRUE, legend = FALSE, 
                                              title = 'Saliva')
pl.plotInd.final.splsda.metab.saliva <- plotInd.final.splsda.metab.saliva$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

final.splsda.metab.saliva.sig <- splsda(df_metab_saliva_sig.X.log.pareto, Y.sub, 
                                    ncomp = optimal.ncomp.splsda.metab.saliva.sig, 
                                    keepX = optimal.keepX.splsda.metab.saliva.sig)


plotIndiv(final.splsda.metab.saliva.sig, comp = c(1,2), 
          group = Y.sub, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, 
          title = ' (a) sPLS-DA on saliva metabolome, comp 1 & 2')


cim(final.splsda.metab.saliva, row.sideColors = Y.sub.color, comp = 1, 
    transpose = TRUE, col.names = TRUE, row.names = FALSE, legend = legend)

cim(final.splsda.metab.saliva.sig, row.sideColors = Y.sub.color, comp = 1, 
    transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.metab.saliva.final <- perf(final.splsda.metab.saliva, 
                                      folds = 10, nrepeat = 50, cpus = 14,
                                      validation = "Mfold", dist = "max.dist", 
                                      progressBar = TRUE)

perf.splsda.metab.saliva.final.sig <- perf(final.splsda.metab.saliva.sig, 
                                       folds = 10, nrepeat = 50, cpus = 14,
                                       validation = "Mfold", dist = "max.dist", 
                                       progressBar = TRUE)


par(mfrow=c(1,3))
plot(perf.splsda.metab.saliva.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.saliva.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.metab.saliva.final$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot(perf.splsda.metab.saliva.final.sig$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.saliva.final.sig$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.metab.saliva.final.sig$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 3', las =2)
par(mfrow=c(1,1))

plotVar(final.splsda.metab.saliva, comp = c(1,2), cex = 3)
plotVar(final.splsda.metab.saliva.sig, comp = c(1,2), cex = 3)

auc.splsda.metab.saliva = auroc(final.splsda.metab.saliva, roc.comp = 1, print = FALSE) 
auc.splsda.metab.saliva.sig = auroc(final.splsda.metab.saliva.sig, roc.comp = 1, print = FALSE) 
splsda.metab.saliva.selectVar.1 <- selectVar(final.splsda.metab.saliva, comp = 1)$name
splsda.metab.saliva.selectVar.2 <- selectVar(final.splsda.metab.saliva, comp = 2)$name
splsda.metab.saliva.selectVar.3 <- selectVar(final.splsda.metab.saliva, comp = 3)$name
#select comp1
df_metab_saliva.X.log.pareto.select <- df_metab_saliva.X.log.pareto[,splsda.metab.saliva.selectVar.1]

#comp1
plotLoadings.final.splsda.metab.saliva <- plotLoadings(final.splsda.metab.saliva, comp = 1, contrib = 'max', method = 'median', title = "saliva")
df.pld.splsda.metab.saliva <- plotLoadings.final.splsda.metab.saliva %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_saliva_identification, by = "ID") 
write.csv(df.pld.splsda.metab.saliva, "metabolome sPLSDA/df.pld.splsda.metab.saliva.csv", row.names = F)
df.pld.splsda.metab.saliva$identification_pld <- ifelse(is.na(df.pld.splsda.metab.saliva$identification), 
                                                       df.pld.splsda.metab.saliva$ID, 
                                                       df.pld.splsda.metab.saliva$identification)

df.pld.splsda.metab.saliva.identified <- df.pld.splsda.metab.saliva %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.saliva.identified$identification_id <- str_trunc(df.pld.splsda.metab.saliva.identified$identification_id, 100)

pld.splsda.metab.saliva.identified <- ggplot(data=df.pld.splsda.metab.saliva.identified, 
                                            aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.splsda.metab.saliva.MS2 <- df.pld.splsda.metab.saliva %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.saliva.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                         df.pld.splsda.metab.saliva.MS2$identification_MS2)
write.csv(df.pld.splsda.metab.saliva.MS2, "metabolome sPLSDA/df.pld.splsda.metab.saliva.MS2.csv", row.names = F)
pld.splsda.metab.saliva.MS2 <- ggplot(data=df.pld.splsda.metab.saliva.MS2, 
                                     aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                         fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.splsda.metab.saliva.single.id <- df.pld.splsda.metab.saliva %>%
  filter(!is.na(identification)) %>%#filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.saliva.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                     df.pld.splsda.metab.saliva.single.id$identification_single_id)
pld.splsda.metab.saliva.single.id <- ggplot(data=df.pld.splsda.metab.saliva.single.id, 
                                           aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                               fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="saliva\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#### serum ####

pca.metab.serum = pca(df_metab_serum.X.log.pareto, ncomp = 10, center = TRUE, scale = T) 
plot(pca.metab.serum)
plotIndiv(pca.metab.serum, group = Y.sub, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on serum metabolome, comp 1 - 2')

pca.metab.serum.sig = pca(df_metab_serum_sig.X.log.pareto, ncomp = 10, center = TRUE, scale = T) 
plot(pca.metab.serum.sig)
plotIndiv(pca.metab.serum.sig, group = Y.sub, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on serum metabolome, comp 1 - 2')

splsda.metab.serum <- splsda(df_metab_serum.X.log.pareto, Y.sub, ncomp = 10)
splsda.metab.serum.sig <- splsda(df_metab_serum_sig.X.log.pareto, Y.sub, ncomp = 10)

plotIndiv(splsda.metab.serum, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')
plotIndiv(splsda.metab.serum.sig, comp = 1:2, 
          group = Y.sub, ind.names = FALSE,  
          ellipse = TRUE, 
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

background.metab.serum = background.predict(splsda.metab.serum, comp.predicted=2, dist = "max.dist")
background.metab.serum.sig = background.predict(splsda.metab.serum.sig, comp.predicted=2, dist = "max.dist")

plotIndiv(splsda.metab.serum, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.metab.serum, 
          legend = TRUE, title = " (b) PLSDA with prediction background")

plotIndiv(splsda.metab.serum.sig, comp = 1:2,
          group = Y.sub, ind.names = FALSE, 
          background = background.metab.serum.sig, 
          legend = TRUE, title = " (b) PLSDA with prediction background")


perf.splsda.metab.serum <- perf(splsda.metab.serum, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

perf.splsda.metab.serum.sig <- perf(splsda.metab.serum.sig, validation = "Mfold", 
                                folds = 10, nrepeat = 50, cpus = 14,
                                progressBar = T, auc = TRUE) 

plot(perf.splsda.metab.serum, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.serum$choice.ncomp

plot(perf.splsda.metab.serum.sig, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.metab.serum.sig$choice.ncomp

list.keepX.splsda.metab.serum <- c(seq(10, 100, 10), seq(150,250,50), seq(300, 500, 100))
list.keepX.splsda.metab.serum.sig <- c(3:9, seq(10, 65, 5))

tune.splsda.metab.serum <- tune.splsda(df_metab_serum.X.log.pareto, Y.sub, ncomp = 3, 
                                       validation = 'Mfold',
                                       folds = 20, nrepeat = 50, 
                                       dist = 'max.dist', 
                                       measure = "BER", 
                                       test.keepX = list.keepX.splsda.metab.serum,
                                       progressBar = T,
                                       cpus = 14) 

tune.splsda.metab.serum.sig <- tune.splsda(df_metab_serum_sig.X.log.pareto, Y.sub, ncomp = 3, 
                                       validation = 'Mfold',
                                       folds = 20, nrepeat = 50, 
                                       dist = 'max.dist', 
                                       measure = "BER", 
                                       test.keepX = list.keepX.splsda.metab.serum.sig,
                                       progressBar = T,
                                       cpus = 14) 

plot(tune.splsda.metab.serum, col = color.jet(3))

tune.splsda.metab.serum$choice.ncomp$ncomp
tune.splsda.metab.serum$choice.keepX

plot(tune.splsda.metab.serum.sig, col = color.jet(3))

tune.splsda.metab.serum.sig$choice.ncomp$ncomp
tune.splsda.metab.serum.sig$choice.keepX

optimal.ncomp.splsda.metab.serum <- tune.splsda.metab.serum$choice.ncomp$ncomp #1
optimal.ncomp.splsda.metab.serum <- 2
optimal.keepX.splsda.metab.serum <- tune.splsda.metab.serum$choice.keepX[1:optimal.ncomp.splsda.metab.serum]

optimal.ncomp.splsda.metab.serum.sig <- tune.splsda.metab.serum.sig$choice.ncomp$ncomp
optimal.keepX.splsda.metab.serum.sig <- tune.splsda.metab.serum.sig$choice.keepX[1:optimal.ncomp.splsda.metab.serum.sig]

final.splsda.metab.serum <- splsda(df_metab_serum.X.log.pareto, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.metab.serum, 
                                   keepX = optimal.keepX.splsda.metab.serum)

final.splsda.metab.serum.sig <- splsda(df_metab_serum_sig.X.log.pareto, Y.sub, 
                                   ncomp = optimal.ncomp.splsda.metab.serum.sig, 
                                   keepX = optimal.keepX.splsda.metab.serum.sig)


plotInd.final.splsda.metab.serum <- plotIndiv(final.splsda.metab.serum, comp = c(1,2), 
                                               group = Y.sub, ind.names = FALSE, 
                                               col.per.group = c("#00468B", "#AD002A"),legend.title = "",
                                               ellipse = TRUE, legend = FALSE, 
                                               title = 'Serum')
pl.plotInd.final.splsda.metab.serum <- plotInd.final.splsda.metab.serum$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

plotIndiv(final.splsda.metab.serum.sig, comp = c(1,2), 
          group = Y.sub, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, 
          title = ' (a) sPLS-DA on serum metabolome, comp 1 & 2')


cim(final.splsda.metab.serum, row.sideColors = Y.sub.color, comp = 1, 
    transpose = TRUE, col.names = TRUE, row.names = FALSE, legend = legend)

cim(final.splsda.metab.serum.sig, row.sideColors = Y.sub.color, comp = 1, 
    transpose = TRUE, col.names = FALSE, row.names = FALSE, legend = legend)

perf.splsda.metab.serum.final <- perf(final.splsda.metab.serum, 
                                       folds = 10, nrepeat = 50, cpus = 14,
                                       validation = "Mfold", dist = "max.dist", 
                                       progressBar = TRUE)

perf.splsda.metab.serum.final.sig <- perf(final.splsda.metab.serum.sig, 
                                           folds = 10, nrepeat = 50, cpus = 14,
                                           validation = "Mfold", dist = "max.dist", 
                                           progressBar = TRUE)


par(mfrow=c(1,2))
plot(perf.splsda.metab.serum.final$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.serum.final$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)

par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(perf.splsda.metab.serum.final.sig$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.metab.serum.final.sig$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
par(mfrow=c(1,1))

plotVar(final.splsda.metab.serum, comp = c(1,2), cex = 3)
plotVar(final.splsda.metab.serum.sig, comp = c(1,2), cex = 3)

auc.splsda.metab.serum = auroc(final.splsda.metab.serum, roc.comp = 1, print = FALSE) 
auc.splsda.metab.serum.sig = auroc(final.splsda.metab.serum.sig, roc.comp = 1, print = FALSE) 
splsda.metab.serum.selectVar <- selectVar(final.splsda.metab.serum, comp = 1)$name

#select comp1
df_metab_serum.X.log.pareto.select <- df_metab_serum.X.log.pareto[,splsda.metab.serum.selectVar]



#comp1
plotLoadings.final.splsda.metab.serum <- plotLoadings(final.splsda.metab.serum, comp = 1, contrib = 'max', method = 'median', title = "serum")
df.pld.splsda.metab.serum <- plotLoadings.final.splsda.metab.serum %>%
  tibble::rownames_to_column("ID") %>% 
  dplyr::left_join(df_metab_serum_identification, by = "ID") 
write.csv(df.pld.splsda.metab.serum, "metabolome sPLSDA/df.pld.splsda.metab.serum.csv", row.names = F)
df.pld.splsda.metab.serum$identification_pld <- ifelse(is.na(df.pld.splsda.metab.serum$identification), 
                                                        df.pld.splsda.metab.serum$ID, 
                                                        df.pld.splsda.metab.serum$identification)

df.pld.splsda.metab.serum.identified <- df.pld.splsda.metab.serum %>%
  filter(!is.na(identification)) %>%#filter out unidentified features by both MS1 and MS2
  mutate(identification_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.serum.identified$identification_id <- str_trunc(df.pld.splsda.metab.serum.identified$identification_id, 100)

pld.splsda.metab.serum.identified <- ggplot(data=df.pld.splsda.metab.serum.identified, 
                                             aes(x=importance, y = reorder(identification_id, -abs(importance)), 
                                                 fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.pld.splsda.metab.serum.MS2 <- df.pld.splsda.metab.serum %>%
  filter(!is.na(MS2Metabolite)) %>%#filter out unidentified features by MS2
  mutate(identification_MS2 = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.serum.MS2$identification_MS2 <- gsub("\\(MS2 NA\\)","(MS2)",
                                                          df.pld.splsda.metab.serum.MS2$identification_MS2)
write.csv(df.pld.splsda.metab.serum.MS2, "metabolome sPLSDA/df.pld.splsda.metab.serum.MS2.csv", row.names = F)
pld.splsda.metab.serum.MS2 <- ggplot(data=df.pld.splsda.metab.serum.MS2, 
                                      aes(x=importance, y = reorder(identification_MS2, -abs(importance)), 
                                          fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.4, 0.4), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df.pld.splsda.metab.serum.single.id <- df.pld.splsda.metab.serum %>%
  filter(!is.na(identification)) %>%#filter out unidentified features 
  filter(!(grepl(';', ident)&is.na(MS2Metabolite))) %>%#filter out multiple identification if no MS2
  mutate(identification_single_id = paste(ID, identification_pld, sep = ": "))
df.pld.splsda.metab.serum.single.id$identification_single_id <- gsub("\\(MS2 NA\\)","(MS2)",
                                                                      df.pld.splsda.metab.serum.single.id$identification_single_id)
pld.splsda.metab.serum.single.id <- ggplot(data=df.pld.splsda.metab.serum.single.id, 
                                            aes(x=importance, y = reorder(identification_single_id, -abs(importance)), 
                                                fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="serum\nmetabolome", fill = "") +
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_x_continuous(limits=c(-0.3, 0.3), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####functional enrichment####
df_metab_feces_Case_Ctrl.all <- read_xlsx("Data filtering/feces/Case_Ctrl.all.xlsx") 
df_metab_saliva_Case_Ctrl.all <- read_xlsx("Data filtering/saliva/Case_Ctrl.all.xlsx") 
df_metab_serum_Case_Ctrl.all <- read_xlsx("Data filtering/serum/Case_Ctrl.all.xlsx") 
df_metab_feces_Case_Ctrl.all.splsda <- df_metab_feces_Case_Ctrl.all %>%
  dplyr::filter(ID %in% row.names(plotLoadings.final.splsda.metab.feces))
df_metab_saliva_Case_Ctrl.all.splsda <- df_metab_saliva_Case_Ctrl.all %>%
  dplyr::filter(ID %in% row.names(plotLoadings.final.splsda.metab.saliva))
df_metab_serum_Case_Ctrl.all.splsda <- df_metab_serum_Case_Ctrl.all %>%
  dplyr::filter(ID %in% row.names(plotLoadings.final.splsda.metab.serum))
df.ttest.feces.all <- read.csv("Data filtering/feces/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)
df.ttest.saliva.all <- read.csv("Data filtering/saliva/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)
df.ttest.serum.all <- read.csv("Data filtering/serum/t_test_all_filt_features.csv") %>%
  rename(ID = X, t.score = t.stat)
#prepare for functional analysis by MetaboAnalyst
df_metab_feces.splsda <- df.ttest.feces.all %>%
  dplyr::right_join(df_metab_feces_Case_Ctrl.all.splsda, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT) %>%
  arrange(mode)
write.csv(df_metab_feces.splsda, "metabolome sPLSDA/df_metab_feces.splsda_for_func.csv", row.names = F)
df_metab_saliva.splsda <- df.ttest.saliva.all %>%
  dplyr::right_join(df_metab_saliva_Case_Ctrl.all.splsda, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT) %>%
  arrange(mode)
write.csv(df_metab_saliva.splsda, "metabolome sPLSDA/df_metab_saliva.splsda_for_func.csv", row.names = F)
df_metab_serum.splsda <- df.ttest.serum.all %>%
  dplyr::right_join(df_metab_serum_Case_Ctrl.all.splsda, by = "ID") %>%
  dplyr::select(ID, MZ, RT ,p.value, t.score) %>% 
  mutate(mode = ifelse(substr(ID, 1, 3) == "pos", "positive", "negative")) %>%
  tibble::column_to_rownames("ID") %>%
  rename(`m.z` = MZ, `r.t` = RT) %>%
  arrange(mode)
write.csv(df_metab_serum.splsda, "metabolome sPLSDA/df_metab_serum.splsda_for_func.csv", row.names = F)

#functional enrichment mummichog and GSEA
library(MetaboAnalystR)
setwd("metabolome sPLSDA")
mSet_fun <-InitDataObjects("mass_all", "mummichog", FALSE)
mSet_fun<-SetPeakFormat(mSet_fun, "rmp")
mSet_fun<-UpdateInstrumentParameters(mSet_fun, 10.0, "mixed", "yes", 0.02)
mSet_fun<-Read.PeakListData(mSet_fun, "df_metab_serum.splsda_for_func.csv")
mSet_fun<-SanityCheckMummichogData(mSet_fun)


mSet_fun<-SetPeakEnrichMethod(mSet_fun, "integ", "v2")
mSet_fun<-SetMummichogPval(mSet_fun, 0.01)
mSet_fun<-PerformPSEA(mSet_fun, "hsa_mfn", "current", 3 , 100)
mSet_fun<-PlotPSEAIntegPaths(mSet_fun, "func_enrich_serum_mummichog2_GSEA", "pdf", 72, width=NA)


####AUC Cliff####
library(ROCR)
library(rcompanion)
library(ggpubr)
#feces
roc.perf.feces <- list()
auc.perf.feces <- list()
df_for_ROC_Cliff<- data.frame(df_metab_feces.X.log.pareto.select, Disease = Y.sub)
df.auc.cliff.feces <- data.frame(matrix(ncol=2,nrow=ncol(df_metab_feces.X.log.pareto.select), 
                                        dimnames=list(colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)]), 
                                                      c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)])){
  pred <- prediction(df_for_ROC_Cliff[,i],df_for_ROC_Cliff$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.feces[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.feces[[i]] <- auc.perf
  df.auc.cliff.feces[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.feces[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                               data = df_for_ROC_Cliff)
}
rownames(df.auc.cliff.feces) <- colnames(df_metab_feces.X.log.pareto.select)
df.auc.cliff.feces.GroupContrib <- merge(df.auc.cliff.feces, 
                                         plotLoadings.final.splsda.metab.feces, by = 0)
df.auc.cliff.feces.GroupContrib$GroupContrib <- factor(df.auc.cliff.feces.GroupContrib$GroupContrib, 
                                                       levels = c("NPD", "PD"))

shapiro.test(df.auc.cliff.feces.GroupContrib$AUC)
shapiro.test(df.auc.cliff.feces.GroupContrib$Cliff)

pl.AUC.feces <- ggboxplot(df.auc.cliff.feces.GroupContrib, 
                          x="GroupContrib", y="AUC", color = "GroupContrib",
                          fill = "GroupContrib", alpha = 0.2,
                          palette = c("#00468B", "#AD002A", "#42B540"),
                          notch = T, add = "jitter") + 
  labs(title = "Feces", x ="Group of Contribution", color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 
pl.Cliff.feces <- ggboxplot(df.auc.cliff.feces.GroupContrib, 
                            x="GroupContrib", y="Cliff", color = "GroupContrib",
                            fill = "GroupContrib", alpha = 0.2,
                            palette = c("#00468B", "#AD002A", "#42B540"),
                            notch = T, add = "jitter") + 
  labs(title = "Feces", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 


#saliva
roc.perf.saliva <- list()
auc.perf.saliva <- list()
df_for_ROC_Cliff<- data.frame(df_metab_saliva.X.log.pareto.select, Disease = Y.sub)
df.auc.cliff.saliva <- data.frame(matrix(ncol=2,nrow=ncol(df_metab_saliva.X.log.pareto.select), 
                                        dimnames=list(colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)]), 
                                                      c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)])){
  pred <- prediction(df_for_ROC_Cliff[,i],df_for_ROC_Cliff$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.saliva[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.saliva[[i]] <- auc.perf
  df.auc.cliff.saliva[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.saliva[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                               data = df_for_ROC_Cliff)
}
rownames(df.auc.cliff.saliva) <- colnames(df_metab_saliva.X.log.pareto.select)
df.auc.cliff.saliva.GroupContrib <- merge(df.auc.cliff.saliva, 
                                         plotLoadings.final.splsda.metab.saliva, by = 0)
df.auc.cliff.saliva.GroupContrib$GroupContrib <- factor(df.auc.cliff.saliva.GroupContrib$GroupContrib, 
                                                       levels = c("NPD", "PD"))

shapiro.test(df.auc.cliff.saliva.GroupContrib$AUC)
shapiro.test(df.auc.cliff.saliva.GroupContrib$Cliff)

pl.AUC.saliva <- ggboxplot(df.auc.cliff.saliva.GroupContrib, 
                          x="GroupContrib", y="AUC", color = "GroupContrib",
                          fill = "GroupContrib", alpha = 0.2,
                          palette = c("#00468B", "#AD002A", "#42B540"),
                          notch = T, add = "jitter") + 
  labs(title = "Saliva", x ="Group of Contribution", color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 
pl.Cliff.saliva <- ggboxplot(df.auc.cliff.saliva.GroupContrib, 
                            x="GroupContrib", y="Cliff", color = "GroupContrib",
                            fill = "GroupContrib", alpha = 0.2,
                            palette = c("#00468B", "#AD002A", "#42B540"),
                            notch = T, add = "jitter") + 
  labs(title = "Saliva", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 



#serum
roc.perf.serum <- list()
auc.perf.serum <- list()
df_for_ROC_Cliff<- data.frame(df_metab_serum.X.log.pareto.select, Disease = Y.sub)
df.auc.cliff.serum <- data.frame(matrix(ncol=2,nrow=ncol(df_metab_serum.X.log.pareto.select), 
                                         dimnames=list(colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)]), 
                                                       c("AUC", "Cliff"))))

for (i in colnames(df_for_ROC_Cliff[1:(ncol(df_for_ROC_Cliff)-1)])){
  pred <- prediction(df_for_ROC_Cliff[,i],df_for_ROC_Cliff$Disease)
  roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
  roc.perf.serum[[i]] <- roc.perf
  auc.perf = performance(pred, measure = "auc")
  auc.perf.serum[[i]] <- auc.perf
  df.auc.cliff.serum[i, "AUC"] <- auc.perf@y.values
  df.auc.cliff.serum[i, "Cliff"] <- cliffDelta(as.formula(paste0(i, " ~ Disease")), 
                                                data = df_for_ROC_Cliff)
}
rownames(df.auc.cliff.serum) <- colnames(df_metab_serum.X.log.pareto.select)
df.auc.cliff.serum.GroupContrib <- merge(df.auc.cliff.serum, 
                                          plotLoadings.final.splsda.metab.serum, by = 0)
df.auc.cliff.serum.GroupContrib$GroupContrib <- factor(df.auc.cliff.serum.GroupContrib$GroupContrib, 
                                                        levels = c("NPD", "PD"))

shapiro.test(df.auc.cliff.serum.GroupContrib$AUC)
shapiro.test(df.auc.cliff.serum.GroupContrib$Cliff)

pl.AUC.serum <- ggboxplot(df.auc.cliff.serum.GroupContrib, 
                           x="GroupContrib", y="AUC", color = "GroupContrib",
                           fill = "GroupContrib", alpha = 0.2,
                           palette = c("#00468B", "#AD002A", "#42B540"),
                           notch = T, add = "jitter") + 
  labs(title = "Serum", x ="Group of Contribution", color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 
pl.Cliff.serum <- ggboxplot(df.auc.cliff.serum.GroupContrib, 
                             x="GroupContrib", y="Cliff", color = "GroupContrib",
                             fill = "GroupContrib", alpha = 0.2,
                             palette = c("#00468B", "#AD002A", "#42B540"),
                             notch = T, add = "jitter") + 
  labs(title = "Serum", x ="Group of Contribution", y= "Cliff's delta", color = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0)) +
  stat_compare_means(label = "p.format", label.y = 0.8) 


lay <- rbind(c(1,2,3),
             c(4,5,6))
legend <- lemon::g_legend(pl.AUC.feces + theme(legend.position='bottom'))
gridExtra::grid.arrange(pl.AUC.saliva + theme(legend.position='hidden'), 
             pl.AUC.serum + theme(legend.position='hidden'), 
             pl.AUC.feces + theme(legend.position='hidden'), 
             pl.Cliff.saliva + theme(legend.position='hidden'), 
             pl.Cliff.serum + theme(legend.position='hidden'),
             pl.Cliff.feces + theme(legend.position='hidden'), 
             bottom=legend$grobs[[1]],
             layout_matrix = lay)

options(ggpubr.parse_aes = FALSE) 

####boxplot of serum sex hormones####

ggboxplot(reshape2::melt(data.frame(df_metab_serum.X.log.pareto[,c("pos-M315T289", 
                                                                   "neg-M367T226", 
                                                                   "neg-M269T185",
                                                                   "neg-M349T222",
                                                                   "neg-M383T194")], 
                     Disease = Y.sub)),
          x = "variable", y = "value",
          color = "Disease", fill = "Disease", alpha = 0.2,
          add = "jitter") +
  labs(x ="", y="log10(Pareto scaling)") + 
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.format") 


####boxplot of significant features####

ggboxplot(reshape2::melt(data.frame(df_metab_feces.X.log.pareto[,c("neg-M593T177", 
                                                                   "pos-M190T161",
                                                                   "neg-M188T159",
                                                                   "neg-M279T329",
                                                                   "pos-M146T46",
                                                                   "neg-M289T160_1",
                                                                   "neg-M287T190",
                                                                   "pos-M289T194")], 
                                    Disease = Y.sub)),
          x = "variable", y = "value",
          color = "Disease", fill = "Disease", alpha = 0.2,
          add = "jitter") +
  labs(x ="", y="log10(Pareto scaling)") + 
  scale_fill_manual(values=c("#00468B", "#AD002A")) +
  scale_color_manual(values=c("#00468B", "#AD002A")) +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) +
  stat_compare_means(aes(group = Disease), method = "wilcox.test", label = "p.format") 



options(ggpubr.parse_aes = FALSE) #let ggpubr to handle non-standard column names



