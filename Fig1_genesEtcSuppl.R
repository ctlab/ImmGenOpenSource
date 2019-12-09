library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Biobase)
source("utils.R")
load("~/Documents/immGen/final_objects/337_es.top12k.Rda")
es.top12k

# genes <- fData(es.top12k)$gene[startsWith(fData(es.top12k)$gene, "Clec")] # 19
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

ggsave('~/Desktop/Fig1_SupplGenes.pdf', width = 16, height = 11, dpi = 300)



# SeqDepth & norm barplots

raw0 <- readr::read_tsv("~/Documents/immGen/OSMNP_unnormalized_genes_count_10_3_18.count_table")
raw <- raw0[, which(colnames(raw0) %in% colnames(es.top12k))]
dim(raw)
all.equal(colnames(raw), colnames(es.top12k))
seq_depth <- apply(raw, 2, sum)

pcaPlot(es.top12k, 1, 2) + aes(color=seq_depth) +
  scale_color_gradientn(colours = gradientColor) +
  ggtitle(paste0("Sequencing depth of the samples \n(based on 10_3_18.count dataset)")) +
  theme_bw()
ggsave('~/Desktop/Fig1_suppl_sDepth.pdf', width = 5, height = 4, dpi = 75) #dpi = 500)



# (preparing two datasets)
# (1) rawRLE
{
  rawRLE <- dplyr::as_tibble(readr::read_csv(
    "/home/octopus/Desktop/immGen/Official_OSMNP_normalized_gene_table.6.19.18.csv"))
  rawRLE <- rawRLE[, -which(colnames(rawRLE) %in% c("MF.6Clo.BM.4",
                                                    "MF.6Clo.BM.5",
                                                    "MF.6Clo.BM.6"))]
  rawRLE$m <- apply(rawRLE, 1, mean)
  rawRLE <- rawRLE[head(order(rawRLE$m, decreasing = T), n=12000), ]
  fdataRLE <- log2(rawRLE[, -c(1,2,ncol(rawRLE))] + 1)
  to_drawRLE <- to_drawRLE[, order(es.top12k$batch)]
  dim(to_drawRLE)
  # all.equal(fdata, to_draw)

  xxx <- cbind(rawRLE[,1], to_drawRLE)
  df_wide <- xxx %>%
    rename(gene = NAME) %>%
    # rename(gene = gene_symbol) %>%
    {
      dups <- duplicated(.$gene) | duplicated(.$gene, fromLast = TRUE)
      if (sum(dups) > 0)
        cat(str_interp('Removed ${sum(dups)} duplicates\n'))
      filter(., !dups)
    } %>%
    select(gene, everything())

  batch_df <- cbind(colnames(df_wide)[-1],
                    es.top12k$batch[order(es.top12k$batch)])

  colnames(batch_df) <- c("sample", "batch")
  batch_df <- as_tibble(batch_df)

  all.equal(batch_df$sample, colnames(df_wide[, -1]))

  df_long <- df_wide %>%
    gather(key = 'sample', value = 'value', -gene)

  # dim(df_long)
  # View(head(df_long))
  df_long <- df_long %>% inner_join(batch_df, by = "sample")
  # View(head(df_long))
  # df_long$batch[565:595]
  # df_long$batch[2037565:2037595]

  # head(df_long)
  # tail(df_long)

  p2 <- df_long %>%
    mutate(sample = fct_reorder(sample, order(batch))) %>%
    ggplot( aes(x=sample, y=value, fill=batch)) +
    scale_fill_brewer(palette="Set1") +
    ggtitle("Log of top 12000 RLE normalized counts") +
    theme(axis.text.x=element_blank()) +
    geom_boxplot()
}


# (2) raw
{
  raw0 <- readr::read_tsv(
    "/home/octopus/Desktop/immGen/OSMNP_unnormalized_genes_count_10_3_18.count_table")

  raw <- raw0[, c(1, which(colnames(raw0) %in% colnames(es.top12k)))]
  dim(raw)

  raw$m <- apply(raw[, -1], 1, mean)
  raw <- raw[head(order(raw$m, decreasing = T), n=12000), ]
  fdata <- log2(raw[, -c(1,ncol(raw))] + 1)
  to_draw <- fdata[, order((es.top12k)$batch)]
  dim(fdata)
  dim(to_draw)
  # all.equal(fdata, to_draw)

  xxx <- cbind(raw[,1], to_draw)
  df_wide <- xxx %>%
    # rename(gene = NAME) %>%
    rename(gene = gene_symbol) %>%
    {
      dups <- duplicated(.$gene) | duplicated(.$gene, fromLast = TRUE)
      if (sum(dups) > 0)
        cat(str_interp('Removed ${sum(dups)} duplicates\n'))
      filter(., !dups)
    } %>%
    select(gene, everything())

  batch_df <- cbind(colnames(df_wide)[-1],
                    es.top12k$batch[order(es.top12k$batch)])

  colnames(batch_df) <- c("sample", "batch")
  batch_df <- as_tibble(batch_df)

  all.equal(batch_df$sample, colnames(df_wide[, -1]))

  df_long <- df_wide %>%
    gather(key = 'sample', value = 'value', -gene)

  # dim(df_long)
  # View(head(df_long))
  df_long <- df_long %>% inner_join(batch_df, by = "sample")
  # View(head(df_long))
  # df_long$batch[565:595]
  # df_long$batch[2037565:2037595]

  # head(df_long)
  # tail(df_long)

  p1 <- df_long %>%
    mutate(sample = fct_reorder(sample, order(batch))) %>%
    ggplot( aes(x=sample, y=value, fill=batch)) +
    scale_fill_brewer(palette="Set1") +
    ggtitle("Log of top 12000 raw counts") +
    theme(axis.text.x=element_blank()) +
    # scale_x_reverse() +
    geom_boxplot()
}



frs_row <- plot_grid(p1, # + theme(legend.position="rigth"),
                     labels = "A",
                     nrow = 1)

seq_row <- plot_grid(p2,
                     labels = "B",
                     nrow = 1)

thd_row <- plot_grid(p3,
                     p4,
                     p5,
                     labels = c("C", "D", "E"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1),
                     rel_heights = c(1, 1, 1),
                     nrow = 1)

q <- plot_grid(frs_row, seq_row, thd_row,
               rel_widths = c(1, 1, 1),
               rel_heights = c(1, 1, 0.7),
               ncol = 1)

ggsave('~/Desktop/immGen/Figs/2_PCpth/Fig2_suppl_qual.pdf', width = 18, height = 18, dpi = 75) #dpi = 500)


