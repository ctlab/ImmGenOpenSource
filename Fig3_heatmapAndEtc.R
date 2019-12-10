session::restore.session("Data/session.RDa")

curRev <- revs[[k]]
expTable <- t(apply(expTable, 1, function (x) mosaic::zscore(x)))
eyegene <- as.data.frame(matrix(nrow=length(curRev$modules), ncol=ncol(gene.exprs)))

for (i in seq_along(curRev$modules)) {
  print(paste("module #", i))
  gs <- unique(E(curRev$modules[[i]])[order(score)]$origin)
  print(length(gs))
  heatmap <- expTable[gs, , drop=F]
  eyegene[i, ] <- apply(heatmap, 2, mean)
}
colnames(eyegene) <- colnames(expTable)
rownames(eyegene) <- paste("Module", 1:9)
View(eyegene)

out <- pheatmap(
  normalize.rows(eyegene),
  cluster_rows=F, cluster_cols=F,
  file=("hierModulesWithColorAnno.pdf"),
  width=27, height=7,
  show_rownames=T, show_colnames=T)

load("Data/337_es.total.Rda")
View(head(pData(es)))
dim(eyegene)
head(eyegene)
annotat <- pData(es)[colnames(eyegene), ]
dim(annotat)
rrr <- data.frame(module = rownames(eyegene))
rownames(rrr) <- rownames(eyegene)
modules <- Biobase::ExpressionSet(as.matrix(eyegene),
                             phenoData = makeAnnotated(annotat),
                             featureData = makeAnnotated(rrr))
ess <- list(modules=modules)
save(ess, file="modules_ohmo8aeLj.rda")
# old one

### PCA on modules space
anno <- rUtils::read.tsv("Data/243_anno.tsv")
head(anno)
View(head(anno))
library(tidyverse)
# *
eyegene2 <- eyegene %>% t %>% as_tibble(rownames = "ImmGen.Nomenclature") %>%
  left_join(anno[, c("ImmGen.Nomenclature", "metaSample")]) %>%
  group_by(metaSample) %>% summarise_each(funs(mean)) %>%
  select(-ImmGen.Nomenclature) %>% t # %>% View()
# *
eyegene2 <- eyegene %>% t %>% as_tibble(rownames = "ImmGen.Nomenclature") %>%
  left_join(anno[, c("ImmGen.Nomenclature", "metaSample")]) %>%
  select(-ImmGen.Nomenclature) %>% t # %>% View()

View(eyegene2)
colnames(eyegene2) <- eyegene2[10, ] # 1 *
eyegene2 <- eyegene2[-10, ] # 1 *
eyegene2 <- apply(eyegene2, 2, as.numeric)
View(eyegene2)
str(eyegene2)
### doing PCA
es <- eyegene2; c1=1; c2=2
pcaPlot(eyegene2, 1, 2, scale = F) + aes(color=colnames(eyegene2)) +
  scale_color_manual(values = metaColors, name = "metaSample") +
  theme_bw()

### doing PCA colored by new clusters
load("Data/337_es.top12k.Rda")
es.top12k$metaSample <- as.character(
  metaSample$value[match(es.top12k$ImmGen.Nomenclature, anno$metaSample)])
es.top12k$metaSample[is.na(es.top12k$metaSample)] <- "treated"
pcaPlot(es.top12k, 1, 2) + aes(color=metaSample) +
  scale_color_manual(values = c(metaColorsList, "grey60"))  + #
  theme_bw()

### dynamic modules representation across all metaSamples
eyegene2 <- eyegene %>% t %>% as_tibble(rownames = "ImmGen.Nomenclature") %>%
  left_join(anno[, c("ImmGen.Nomenclature", "metaSample")]) %>%
  group_by(metaSample) %>% summarise_each(funs(mean)) %>%
  dplyr::select(-ImmGen.Nomenclature) # %>% View()
eyegene2 <- eyegene2[c(2, 5, 4, 10, 3, 6, 7, 8, 9, 1), ]

for(i in seq_along(1:10)){
  print(i)
  data <- eyegene2[, c(1, i + 1)]

  data[, 2] <- scales::rescale(unlist(unname(eyegene2[, i + 1])), to=c(0,1))

  unname(unlist(eyegene2[, 1]))
  nnn <- c("EMP (E10.5)", "MF (E6-E8)", "alv MF & SPM", "tissue MF", "adipose MF", "MG", "Mo", "pDC", "tissue DC", "mig DC")
  data[, 1] <- factor(nnn, levels = nnn)

  ggplot(data = data, aes(x = as.factor(data$metaSample), y = unlist(data[, 2]), group=1)) +
    geom_line(color = "grey60", size = 2) + geom_point(size = 3) +
    theme(axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust=1))
}

