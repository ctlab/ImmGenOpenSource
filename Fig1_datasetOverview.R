source("utils.R")
library(plotly)
load("Data/337_es.top12k.Rda")

# es.top12k$Cell.Family[es.top12k$Organ == "Yolk sac"] <- "EMP"
# es.top12k$Cell.Family[grep("pDC", es.top12k$ImmGen.Nomenclature)] <- "pDC"
# es.top12k$Cell.Family[grep("cDC1", es.top12k$ImmGen.Nomenclature)] <- "cDC1"
# es.top12k$Cell.Family[grep("cDC2", es.top12k$ImmGen.Nomenclature)] <- "cDC2"
# es.top12k$Cell.Family[grep("E1[468]", es.top12k$ImmGen.Nomenclature)] <- "MF (E14.5-E18.5)"
# es.top12k$Cell.Family[grep("EB", es.top12k$ImmGen.Nomenclature)] <- "MF (E6-E8)"
# es.top12k$Cell.Family[es.top12k$Cell.Family == "Macrophage"] <- "MF"
# es.top12k$Cell.Family[es.top12k$Cell.Family == "Dendritic Cell"] <- "DC"
# es.top12k$Cell.Family[es.top12k$Cell.Family == "Monocyte"] <- "Mo"
# es.top12k$Cell.Family[es.top12k$Cell.Family == "Microglia"] <- "MG"
# OR
es.top12k$Cell.Family[es.top12k$Organ == "Yolk sac"] <- "Erythro-myeloid progenitor"
es.top12k$Cell.Family[grep("EB", es.top12k$ImmGen.Nomenclature)] <- "Macrophage (E6-E8)"

pcaPlot(es.top12k, 1, 2) + aes(color=es.top12k$Cell.Family) +
  scale_color_manual(values = cellColorsList) +
  theme_bw()

pcaPlot(es.top12k, 1, 2) + aes(color=es.top12k$Organ) +
  scale_color_manual(values = tissueColorsList) +
  theme_bw()

pcaPlot(es.top12k, 1, 2) + aes(color=as.factor(es.top12k$Batch)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set1")[-6], "black")) +
  theme_bw()

pcaPlot(es.top12k, 1, 2) + aes(color=es.top12k$Lab) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(12, "Paired")[-11], "black", "turquoise")) +
  theme_bw()

es.top12k$treatment_var <- "non-stimulated"
es.top12k$treatment_var[es.top12k$treatment != "none"] <- "stimulated"
pcaPlot(es.top12k, 1, 2) + aes(color=es.top12k$treatment_var) +
  scale_color_manual(values = c("grey20", "grey60")) +
  theme_bw()

### make interactive content
pca <- prcomp(t(exprs(es.top12k)))
explained <- (pca$sdev)^2/sum(pca$sdev^2)
xs <- sprintf("PC%s", seq_along(explained))
xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)
a <- cbind(as.data.frame(pca$x), the_sample = colnames(es.top12k),
           pData(es.top12k))
pp <- plot_ly(a, x = ~a$PC1, y = ~a$PC2, type = 'scatter', mode = 'markers',
              color = ~pData(es.top12k)$Cell.Family,
              colors = unlist(unname(cellColorsList))[c(1,2,5,3,4,6)],
              marker = list(size=10),
              text = ~paste0('Sample: ', pData(es.top12k)$ImmGen.Nomenclature,
                            '\nOrgan: ', pData(es.top12k)$Organ,
                            '\nSortMark: ', pData(es.top12k)$Sorting.Markers,
                            '\nStimulation: ', pData(es.top12k)$treatment,
                            '\nLab: ', pData(es.top12k)$Lab,
                            '\nOwner: ', pData(es.top12k)$Owner
                            # '\nBatch: ', pData(es.top12k)$Batch
                            ))
pp
htmlwidgets::saveWidget(pp, file = 'SupplData/ESuppl_datasetOverview.html')




