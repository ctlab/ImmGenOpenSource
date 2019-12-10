source("utils.R")
load("Data/337_es.top12k.Rda")

genes <- c(
  "Lyz2", "H2-Aa", "Ccr2", "Cx3cr1",
  "Mertk", "Lyve1", "Sall2", "Siglecf", # "Adgre1", "Itgax", "Itgam",
  "Zbtb46", "Xcr1", "Sirpa", "Ccr7")
length(genes)

fff <- list()
for(geneNum in seq_along(genes)){
  geneExpres <- exprs(es.top12k)[which(fData(es.top12k)$gene == genes[geneNum]), ]
  fff[[geneNum]] <- scales::rescale(geneExpres, to = c(-2,2))
}

{
  p1 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[1]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[1])) +
    theme_bw()

  p2 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[2]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[2])) +
  theme_bw()

  p3 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[3]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[3])) +
  theme_bw()

  p4 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[4]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[4])) +
  theme_bw()

  p5 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[5]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[5])) +
  theme_bw()

  p6 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[6]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[6])) +
  theme_bw()

  p7 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[7]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[7])) +
  theme_bw()

  p8 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[8]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[8])) +
  theme_bw()

  p9 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[9]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[9])) +
  theme_bw()

  p10 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[10]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[10])) +
  theme_bw()

  p11 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[11]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[11])) +
  theme_bw()

  p12 <- pcaPlot(es.top12k, 1, 2) + aes(color=fff[[12]]) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle(paste0("Expression of ", genes[12])) +
  theme_bw()
}

frs_row <- plot_grid(p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     labels = c("A", "B", "C", "D"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1, 1), nrow = 1)

sec_row <- plot_grid(p5 + theme(legend.position="none"),
                     p6 + theme(legend.position="none"),
                     p7 + theme(legend.position="none"),
                     p8 + theme(legend.position="none"),
                     labels = c("E", "F", "G", "H"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1, 1), nrow = 1)

thd_row <- plot_grid(p9 + theme(legend.position="none"),
                     p10 + theme(legend.position="none"),
                     p11 + theme(legend.position="none"),
                     p12 + theme(legend.position="none"),
                     labels = c("I", "J", "K", "L"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1, 1), nrow = 1)

q <- plot_grid(frs_row, sec_row, thd_row,
               rel_widths = c(1, 1, 1),
               rel_heights = c(1, 1, 1),
               ncol = 1)

###
raw0 <- readr::read_tsv("Data/OSMNP_unnormalized_genes_count_10_3_18.count_table")
raw <- raw0[, which(colnames(raw0) %in% colnames(es.top12k))]
dim(raw)
all.equal(colnames(raw), colnames(es.top12k))
seq_depth <- apply(raw, 2, sum)

pcaPlot(es.top12k, 1, 2) + aes(color=seq_depth) +
  scale_color_gradientn(colours = gradientColor) +
  theme_bw()
