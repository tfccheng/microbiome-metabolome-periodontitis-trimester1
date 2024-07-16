library(Boruta)
library(randomForest)
library(caret)
library(MLeval)
library(phyloseq)
library(mia)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library("readxl")
setwd("")
set.seed(123)
tse.micro <- readRDS("")
sample_meta <- data.frame(colData(tse.micro))


tse.micro.subsamples <- list()
for (n in c("SBP", "MSA", "MST")){
  tse.micro.subsample <- tse.micro[, tse.micro$Sample %in% n] %>%
    subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05)  #filter with 0.05% detection & 5% prevalence
  for (m in c("Genus", "Species")){
    altExp(tse.micro.subsample, m) <- agglomerateByRank(tse.micro.subsample, rank = m, onRankOnly = T)
  }
  tse.micro.subsample <- tse.micro.subsample %>%
    transformSamples(assay_name = "counts", method = "relabundance") %>%
    transformSamples(assay_name = "relabundance", method = "rclr")
  tse.micro.subsamples[[n]] <- tse.micro.subsample
}


########genus level
#### Get assay rclr####
df.SBP.g.rclr <- altExp(tse.micro.subsamples[["SBP"]], "Genus") %>% 
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>%
  assay("rclr") %>% t()%>% as.data.frame()
g.SBP.rclr.names <- colnames(df.SBP.g.rclr)
df.SBP.g.rclr$Disease <- as.factor(sample_meta.SBP$Disease)
names(df.SBP.g.rclr) <- make.names(names(df.SBP.g.rclr))

df.MSA.g.rclr <- altExp(tse.micro.subsamples[["MSA"]], "Genus") %>% 
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>%
  assay("rclr") %>% t()%>% as.data.frame()
g.MSA.rclr.names <- colnames(df.MSA.g.rclr)
df.MSA.g.rclr$Disease <- as.factor(sample_meta.MSA$Disease)
names(df.MSA.g.rclr) <- make.names(names(df.MSA.g.rclr))

df.MST.g.rclr <- altExp(tse.micro.subsamples[["MST"]], "Genus") %>% 
  transformSamples(assay_name = "counts", method = "relabundance") %>%
  transformSamples(assay_name = "relabundance", method = "rclr") %>%
  assay("rclr") %>% t()%>% as.data.frame()
g.MST.rclr.names <- colnames(df.MST.g.rclr)
df.MST.g.rclr$Disease <- as.factor(sample_meta.MST$Disease)
names(df.MST.g.rclr) <- make.names(names(df.MST.g.rclr))


#####SBP caret and Boruta####
set.seed(123)

#10 folds repeat 5 times
rf.control <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 5, allowParallel = T,
                           search='grid')

#Create training grid for mtry values.
rf.tunegrid <- expand.grid(.mtry = 2:sqrt(ncol(df.SBP.g.rclr)))

#Run RandomForests. 

rf.gridsearch.SBP.g.rclr <- train(Disease ~ ., data=df.SBP.g.rclr, method = "rf",     
                             metric = 'Accuracy',
                             trControl= rf.control,
                             tuneGrid = rf.tunegrid)
print(rf.gridsearch.SBP.g.rclr)
plot(rf.gridsearch.SBP.g.rclr)


rf.control <- trainControl(method="repeatedcv", 
                           number=10, 
                           repeats=5, allowParallel = T,
                           search="grid")

rf.tunegrid <- expand.grid(.mtry=6)
modellist.SBP.g.rclr <- list()
for (ntree in c(500, 1000,2000, 3000, 4000, 5000)) {
  set.seed(123) 
  fit <- train(Disease ~ ., data=df.SBP.g.rclr, method="rf", 
               tuneGrid=rf.tunegrid, trControl=rf.control, 
               ntree=ntree)
  key <- toString(ntree)
  modellist.SBP.g.rclr[[key]] <- fit
}

results <- resamples(modellist.SBP.g.rclr)
dotplot(results) 
summary(results)
rf.model.varImp <- varImp(modellist.SBP.g.rclr[["500"]])
imps <- rf.model.varImp$importance
ggplot(varImp(modellist.SBP.g.rclr[["500"]]), top = 20)
plot(varImp(modellist.SBP.g.rclr[["500"]]), top = 20)

df.SBP.g.rclr.40 <- df.SBP.g.rclr %>%
  dplyr::select(-rownames(imps %>% filter(imps$Overall < 40)))
write.csv(imps %>% filter(imps$Overall > 40), "Machine learning/df.SBP.g.rclr.40.csv", row.names = T)
set.seed(123)
boruta_SBP.g.rclr<-Boruta(Disease ~ ., data=df.SBP.g.rclr.40, ntree =500, mtry=6,
          mcAdj = TRUE, doTrace = 2,maxRuns = 500,
          pValue=0.05)
print(boruta_SBP.g.rclr)
par(mar = c(10, 5, 0.1, 0.5))
plot(boruta_SBP.g.rclr,xlab="",las=2,
     cex.axis = 0.7, colCode = c("#00A087", "#EFC000", "#E64B35", "#999999"))
bor <- TentativeRoughFix(boruta_SBP.g.rclr)
print(bor)


####
#####MSA caret and Boruta####
set.seed(123)

#10 folds repeat 3 times
rf.control <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 5, allowParallel = T,
                           search='grid')

#Create training grid for mtry values.
rf.tunegrid <- expand.grid(.mtry = 2:sqrt(ncol(df.MSA.g.rclr)))

#Run RandomForests. 

rf.gridsearch.MSA.g.rclr <- train(Disease ~ ., data=df.MSA.g.rclr, method = "rf",     
                             metric = 'Accuracy',
                             trControl= rf.control,
                             tuneGrid = rf.tunegrid)
print(rf.gridsearch.MSA.g.rclr)
plot(rf.gridsearch.MSA.g.rclr)


rf.control <- trainControl(method="repeatedcv", 
                           number=10, 
                           repeats=5, allowParallel = T,
                           search="grid")

rf.tunegrid <- expand.grid(.mtry=5)
modellist.MSA.g.rclr <- list()
for (ntree in c(500, 1000,2000, 3000, 4000, 5000)) {
  set.seed(123) 
  fit <- train(Disease ~ ., data=df.MSA.g.rclr, method="rf", 
               tuneGrid=rf.tunegrid, trControl=rf.control, 
               ntree=ntree)
  key <- toString(ntree)
  modellist.MSA.g.rclr[[key]] <- fit
}

results <- resamples(modellist.MSA.g.rclr)
dotplot(results) 
summary(results)
rf.model.varImp <- varImp(modellist.MSA.g.rclr[["1000"]])
imps <- rf.model.varImp$importance
ggplot(varImp(modellist.MSA.g.rclr[["1000"]]), top = 20)
plot(varImp(modellist.MSA.g.rclr[["1000"]]), top = 20)


df.MSA.g.rclr.40 <- df.MSA.g.rclr %>%
  dplyr::select(-rownames(imps %>% filter(imps$Overall < 40)))
write.csv(imps %>% filter(imps$Overall > 40), "Machine learning/df.MSA.g.rclr.40.csv", row.names = T)
set.seed(123)
boruta_MSA.g.rclr<-Boruta(Disease ~ ., data=df.MSA.g.rclr.40, ntree =1000, mtry=5,
                     mcAdj = TRUE, doTrace = 2,maxRuns = 500,
                     pValue=0.05)
print(boruta_MSA.g.rclr)
par(mar = c(10, 5, 0.1, 0.5))
plot(boruta_MSA.g.rclr,xlab="",las=2,
     cex.axis = 0.7, colCode = c("#00A087", "#EFC000", "#E64B35", "#999999"))
bor <- TentativeRoughFix(boruta_MSA.g.rclr)
print(bor)


####
#####MST caret and Boruta####
set.seed(123)

#10 folds repeat 3 times
rf.control <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 5, allowParallel = T,
                           search='grid')

#Create training grid for mtry values.
rf.tunegrid <- expand.grid(.mtry = 2:sqrt(ncol(df.MST.g.rclr)))

#Run RandomForests. 

rf.gridsearch.MST.g.rclr <- train(Disease ~ ., data=df.MST.g.rclr, method = "rf",     
                             metric = 'Accuracy',
                             trControl= rf.control,
                             tuneGrid = rf.tunegrid)
print(rf.gridsearch.MST.g.rclr)
plot(rf.gridsearch.MST.g.rclr)


rf.control <- trainControl(method="repeatedcv", 
                           number=10, 
                           repeats=5, allowParallel = T,
                           search="grid")

rf.tunegrid <- expand.grid(.mtry=8)
modellist.MST.g.rclr <- list()
for (ntree in c(500, 1000,2000, 3000, 4000, 5000)) {
  set.seed(123) 
  fit <- train(Disease ~ ., data=df.MST.g.rclr, method="rf", 
               tuneGrid=rf.tunegrid, trControl=rf.control, 
               ntree=ntree)
  key <- toString(ntree)
  modellist.MST.g.rclr[[key]] <- fit
}

results <- resamples(modellist.MST.g.rclr)
dotplot(results) 
summary(results)
rf.model.varImp <- varImp(modellist.MST.g.rclr[["500"]])
imps <- rf.model.varImp$importance
ggplot(varImp(modellist.MST.g.rclr[["500"]]), top = 20)
plot(varImp(modellist.MST.g.rclr[["500"]]), top = 20)

df.MST.g.rclr.40 <- df.MST.g.rclr %>%
  dplyr::select(-rownames(imps %>% filter(imps$Overall < 40)))
# df.MST.g.rclr.50 <- df.MST.g.rclr %>%
#   dplyr::select(-rownames(imps %>% filter(imps$Overall < 50)))
write.csv(imps %>% filter(imps$Overall > 40), "Machine learning/df.MST.g.rclr.40.csv", row.names = T)
set.seed(123)
boruta_MST.g.rclr<-Boruta(Disease ~ ., data=df.MST.g.rclr.40, ntree =500, mtry=8,
                     mcAdj = TRUE, doTrace = 2,maxRuns = 500,
                     pValue=0.05)
print(boruta_MST.g.rclr)
par(mar = c(10, 5, 0.1, 0.5))
plot(boruta_MST.g.rclr,xlab="",las=2,
     cex.axis = 0.7, colCode = c("#00A087", "#EFC000", "#E64B35", "#999999"))
bor <- TentativeRoughFix(boruta_MST.g.rclr)
print(bor)


