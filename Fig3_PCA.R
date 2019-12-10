source("utils.R")
m_dir <-"Data"
m_files <- list.files(m_dir,"m\\.[0-9]+\\.genes\\.tsv", full.names = T)
m_files <- m_files[-2]
load("Data/337_es.top12k.Rda")

zFdata <- t(apply(exprs(es.top12k), 1, function (x) mosaic::zscore(x)))
moduleColors <- list()

for (i in 1:length(m_files)){
  file <- data.table::fread(m_files[i])
  name <- basename(m_files[i])
  name <- sapply(strsplit(name, ".", fixed = T),"[[", 2)
  moduleGenesEntrez <- file$Entrez
  print(length(moduleGenesEntrez))

  fdata <- zFdata[which(fData(es.top12k)$entrez %in% moduleGenesEntrez), ]

  moduleGenes <- apply(fdata, 2, mean)
  moduleGenes[which(pData(es.top12k)$treatment != "none")] <- NA
  moduleColors[[i]] <- moduleGenes
}

{
  p1 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[1]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 1") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p2 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[2]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 2") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p3 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[3]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 4") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p4 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[4]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 5") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p5 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[5]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 6") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p6 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[6]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 7") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p7 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[7]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 3") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p8 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[8]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 8") +
    theme_bw() + theme(plot.title = element_text(size = 25))
  p9 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(moduleColors[[9]], to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Module 9") +
    theme_bw() + theme(plot.title = element_text(size = 25)) #, face = "bold"))
}

f_row <- plot_grid(p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p3 + theme(legend.position="none"),
                   # labels = c("A", "B", "C"),
                   labels = c("", "", ""),
                   hjust = -2,
                   rel_widths = c(1, 1, 1), nrow = 1)

s_row <- plot_grid(p4 + theme(legend.position="none"),
                   p5 + theme(legend.position="none"),
                   p6 + theme(legend.position="none"),
                   # labels = c("D", "E", "F"),
                   labels = c("", "", ""),
                   hjust = -2,
                   rel_widths = c(1, 1, 1), nrow = 1)

t_row <- plot_grid(p7 + theme(legend.position="none"),
                   p8 + theme(legend.position="none"),
                   p9 + theme(legend.position="none"),
                   # labels = c("G", "H", "I"),
                   labels = c("", "", ""),
                   hjust = -2,
                   rel_widths = c(1, 1, 1), nrow = 1)

q <- plot_grid(f_row, s_row, t_row, rel_heights = c(1, 1, 1), ncol = 1)