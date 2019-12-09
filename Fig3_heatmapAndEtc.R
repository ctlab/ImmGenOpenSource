session::restore.session("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/session.RDa")

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



# making color annotation
load("~/Documents/immGen/final_objects/337_es.top12k.Rda")
anno <- pData(es.top12k)[, c("Cell.Family", "Organ")]
# * ordering  like in phantasus
# load("~/Documents/immGen/phantasus_files/ImmGen_total_Eduw0mei4.rda")
load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig56/backup/ImmGen_norm.top12k_dj8dUkl1.rda")
order <- readLines("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig56/samples_order.txt")
ordDat <- data.frame(order = order, pos = 1:length(order))
View(ordDat)
pData(ess$immGen_norm.top12k)$hierSort <- ordDat$pos[match(pData(ess$immGen_norm.top12k)$sample, ordDat$order)]
pData(ess$immGen_norm.top12k)$hierSort[is.na(pData(ess$immGen_norm.top12k)$hierSort)] <- 0
ess <- list(immGen_norm.top12k=ess$immGen_norm.top12k)
save(ess, file="~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig56/backup/ImmGen_norm.top12k_dj8dUkl1.rda")
rUtils::write.gct(ess$immGen_norm.top12k, file="~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig56/backup/1g.gct")
# *
gct <- read.table("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig56/samples_order.txt"); gct$V1
anno <- anno[match(gct$V1, rownames(anno)), ]
eyegene <- eyegene[, match(gct$V1, colnames(eyegene))]

annotation_colors <- list(Organ = unlist(tissueColorsList),
                          Cell.Family = unlist(cellColorsList))

out <- pheatmap(
  normalize.rows(eyegene),
  cluster_rows=F, cluster_cols=F,
  file=("~/Desktop/hierModulesWithColorAnno.pdf"),
  width=27, height=7,
  show_rownames=T, show_colnames=T,
  # clustering_distance_cols = "manhattan",
  # clustering_distance_cols = "euclidean",
  # clustering_method = "average",
  # cutree_cols = 11,
  annotation = anno,
  annotation_colors = annotation_colors)



# doing supplementary
# spyderHolders <- sort(cutree(out$tree_col, k=11))
# write.table(spyderHolders,
#             "~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/metaSamples_12.xlsx",
#             sep="\t")
# save(spyderHolders,
#      file="/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/spyderHolders_12.rda")


# saving modules for phantasus
load("~/Documents/immGen/final_objects/337_es.total_Ph.Rda")
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
save(ess, file="/home/octopus/Documents/immGen/phantasus_files/modules_ohmo8aeLj.rda")
# old one
# save(ess, file="/home/octopus/R-studio/nclust/immGene/new3/ImmGen_10norm.metaSamples.rda")



# PCA on modules space
anno <- rUtils::read.tsv(
  "~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/243_anno.tsv")
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
ggsave("~/Desktop/PCA_on_9_modules_scaleF.pdf", width = 7, height = 4) # width = 7, height = 3)



# doing PCA colored by new clusters
metaSample <- read_csv("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/anno.csv")
load("~/Documents/immGen/final_objects/337_es.top12k.Rda")
es.top12k$metaSample <- as.character(
  metaSample$value[match(es.top12k$ImmGen.Nomenclature, metaSample$samples)])
es.top12k$metaSample[is.na(es.top12k$metaSample)] <- "treated"
# save(es.top12k,
#      file = "~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/es.top12k.Rda")
pcaPlot(es.top12k, 1, 2) + aes(color=metaSample) +
  scale_color_manual(values = c(metaColorsList, "grey60"))  + #
  theme_bw()
ggsave('~/Desktop/PCA_metaSample.pdf', width = 7.5, height = 4.5, dpi = 600)




# prepare data for radar charts drawing
library(tidyverse)
a <- tibble::as.tibble(colnames(expTable))
colnames(a) <- "samples"
b <- tibble::as.tibble(spyderHolders) %>% mutate(samples = names(spyderHolders))
View(a)
a <- a %>% full_join(b, by = "samples")

a$value[is.na(a$value)] <- "stimulated"

rUtils::write.tsv(a,
                  file = "/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/metaSamples_12.tsv")

View(eyegene)
apply(eyegene, 1, range)
apply(normalize.rows(eyegene), 1, range)
eyegene <- normalize.rows(eyegene)
View(a)
save(a,
     file = "/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/a.rda")
save(eyegene,
     file = "/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/eyegene.rda")




# dynamic modules representation across all metaSamples
anno <- rUtils::read.tsv(
  "~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/oneMore_anno.tsv")
head(anno)
View(head(anno))
library(tidyverse)
# *
eyegene2 <- eyegene %>% t %>% as_tibble(rownames = "ImmGen.Nomenclature") %>%
  left_join(anno[, c("ImmGen.Nomenclature", "metaSample")]) %>%
  group_by(metaSample) %>% summarise_each(funs(mean)) %>%
  dplyr::select(-ImmGen.Nomenclature) # %>% View()
View(eyegene2)
eyegene2 <- eyegene2[c(2, 5, 4, 10, 3, 6, 7, 8, 9, 1), ]
str(eyegene2)
View(eyegene2)

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

  ggsave(paste0("~/Desktop/", i, "module.pdf"), width = 6, height = 2.5)
}

