library(ggpubr)
library(paletteer)
setwd("")
feces_MS2superclass <- readRDS("feces_MS2superclass.rds") %>% as.vector.data.frame()
saliva_MS2superclass <- readRDS("saliva_MS2superclass.rds") %>% as.vector.data.frame()
serum_MS2superclass <- readRDS("serum_MS2superclass.rds") %>% as.vector.data.frame()
df <- qpcR::cbind.na(feces = feces_MS2superclass$MS2superclass, 
                      saliva = saliva_MS2superclass$MS2superclass, 
                      serum = serum_MS2superclass$MS2superclass) %>%
  as.data.frame()


feces.count <- data.frame(unclass(summary(factor((df %>% filter(feces != "NA"))$feces)))) %>%
  tibble::rownames_to_column("group")
colnames(feces.count) <- c("group", "feces.count")
saliva.count <- data.frame(unclass(summary(factor((df%>% filter(!is.na(saliva)))$saliva)))) %>%
  tibble::rownames_to_column("group")
colnames(saliva.count) <- c("group", "saliva.count")
serum.count <- data.frame(unclass(summary(factor((df%>% filter(!is.na(serum)))$serum)))) %>%
  tibble::rownames_to_column("group")
colnames(serum.count) <- c("group", "serum.count")

df.count <- dplyr::full_join(feces.count, saliva.count, by = "group")
df.count <- dplyr::full_join(df.count, serum.count, by = "group") %>%
  replace(is.na(.), 0)
df.count$feces.percent <- df.count$feces.count/sum(df.count$feces.count)
df.count$saliva.percent <- df.count$saliva.count/sum(df.count$saliva.count)
df.count$serum.percent <- df.count$serum.count/sum(df.count$serum.count)

feces.limma.sig <- read.csv("df.covariate.BMI_w.anno.feces.csv") %>%
  filter(MS2superclass != "Unknown") %>% dplyr::select(MS2superclass)
feces.limma.sig.count <- data.frame(unclass(summary(factor(feces.limma.sig$MS2superclass)))) %>%
  tibble::rownames_to_column("group")
colnames(feces.limma.sig.count) <- c("group", "feces.count")
saliva.limma.sig <- read.csv("df.covariate.BMI_w.anno.saliva.csv") %>%
  filter(MS2superclass != "Unknown") %>% dplyr::select(MS2superclass)
saliva.limma.sig.count <- data.frame(unclass(summary(factor(saliva.limma.sig$MS2superclass)))) %>%
  tibble::rownames_to_column("group")
colnames(saliva.limma.sig.count) <- c("group", "saliva.count")
serum.limma.sig <- read.csv("df.covariate.BMI_w.anno.serum.csv") %>%
  filter(MS2superclass != "Unknown") %>% dplyr::select(MS2superclass)
serum.limma.sig.count <- data.frame(unclass(summary(factor(serum.limma.sig$MS2superclass)))) %>%
  tibble::rownames_to_column("group")
colnames(serum.limma.sig.count) <- c("group", "serum.count")
df.limma.sig.count <- dplyr::full_join(feces.limma.sig.count, saliva.limma.sig.count, by = "group")
df.limma.sig.count <- dplyr::full_join(df.limma.sig.count, serum.limma.sig.count, by = "group") %>%
  replace(is.na(.), 0)
df.limma.sig.count$feces.percent <- df.limma.sig.count$feces.count/sum(df.limma.sig.count$feces.count)
df.limma.sig.count$saliva.percent <- df.limma.sig.count$saliva.count/sum(df.limma.sig.count$saliva.count)
df.limma.sig.count$serum.percent <- df.limma.sig.count$serum.count/sum(df.limma.sig.count$serum.count)


ggdonutchart(df.count, "feces.count", label = NULL, 
             fill = "group", color = "white", palette = paletteer_d("ggthemes::Classic_20")) 
ggdonutchart(df.count, "saliva.count", label = NULL, 
             fill = "group", color = "white", palette = paletteer_d("ggthemes::Classic_20")) 
ggdonutchart(df.count, "serum.count", label = NULL, 
             fill = "group", color = "white", palette = paletteer_d("ggthemes::Classic_20")) 

df1 <- melt(df.count[, c(1, 7,6,5)])

.remove_axis <- function(){
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    # panel.border = element_blank(),
    # panel.grid=element_blank(),
    axis.ticks = element_blank()
  )
}

ggplot(df1, aes(x = variable, y = value, fill = group))+
  geom_bar(stat = "identity", color = "white", linewidth = 1)+
  #geom_text(aes(label = variable)) +
  coord_polar(theta = "y", start = 0, clip = "off") +
  theme_pubr() + .remove_axis() +
  scale_fill_paletteer_d("ggthemes::Classic_20")

df2 <- melt(df.limma.sig.count[, c(1, 7,6,5)])
ggplot(df2, aes(x = variable, y = value, fill = group))+
  geom_bar(stat = "identity", color = "white", linewidth = 1)+
  #geom_text(aes(label = variable)) +
  coord_polar(theta = "y", start = 0, clip = "off") +
  theme_pubr() + .remove_axis() +
  scale_fill_manual(values = paletteer_d("ggthemes::Classic_20")[c(1,2,5:11,13,18)])

