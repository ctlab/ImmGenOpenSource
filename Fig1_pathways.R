source("utils.R")
load("Data/pathway2name.rda")
load("Data/pathways.rda")
load("Data/337_es.total.Rda")
load("Data/337_es.top12k.Rda")
es
es.top12k

fData(es)$Entrez <- fData(es)$entrez
exprs(es) <- t(apply(exprs(es), 1, function (x) mosaic::zscore(x)))
fdata <- exprs(es)

# mitoNumbers <- grep("^mt", fData(es)$gene)
# fData(es)$gene[mitoNumbers]
# mito <- apply(fdata[mitoNumbers, ], 2, mean)
#
# pcaPlot(es.top12k, 1, 2) + aes(color=scale(mito)) +
#   scale_color_gradientn(colours = gradientColor) + theme_bw()

# pathway2name$PATHNAME[grep("folate", pathway2name$PATHNAME)]
{
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Glycolysis")]
  gP <- pathways$`R-MMU-70171`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Pentose phosphate pathway")]
  pppP <- pathways$`mmu00030`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Citric acid cycle (TCA cycle)")]
  tcaP <- pathways$`R-MMU-71403`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Mitochondrial Fatty Acid Beta-Oxidation")]
  betaP <- pathways$`R-MMU-77289`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Valine, leucine and isoleucine degradation")]
  bcaaP <- pathways$`mmu00280`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Cholesterol biosynthesis")]
  chP <- pathways$`R-MMU-191273`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Biosynthesis of amino acids")]
  aabP <- pathways$`mmu01230`; aabP
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Folate biosynthesis")]
  fP <- pathways$`mmu00790`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Metabolism of folate and pterines")]
  fpP <- pathways$`R-MMU-196757`

  pathway2name$PATHID[which(pathway2name$PATHNAME == "Signaling by Interleukins")]
  ilP <- pathways$`R-MMU-449147`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "Antigen processing-Cross presentation")]
  mhcP <- pathways$`R-MMU-1236975`
  pathway2name$PATHID[which(pathway2name$PATHNAME == "NF-kappa B signaling pathway")]
  nfkbP <- pathways$`mmu04064`
}

{
  # pathway2name$PATHNAME[pathway2name$PATHID == "mmu00790"]
  # length(pathways$`mmu00790`)
  # length(which(fData(es)$Entrez %in% pathways$`mmu00790`))
  gF <- apply(fdata[which(fData(es)$Entrez %in% pathways$`mmu00790`), ], 2, mean)
  # title = "Glycolysis"
  gP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-70171`), ], 2, mean)
  # title = "TCA"
  tcaP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-71403`), ], 2, mean)
  # title = "Mitochondrial Fatty Acid Beta-Oxidation"
  betaP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-77289`), ], 2, mean)
  # title = "VAL, LEU and ILE degradation"
  bcaaP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`mmu00280`), ], 2, mean)
  # title = "Cholesterol biosynthesis"
  chP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-191273`), ], 2, mean)
  # title = "Biosynthesis of amino acids"
  aabP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`mmu01230`), ], 2, mean)
  # etcP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-1428517`), ], 2, mean)
  # title = "Adherens junctions interactions"
  # adhP <- apply(fdata[which(fData(es)$Entrez %in% pathways$`R-MMU-418990`), ], 2, mean)
  # Antigen processing - cross presentation
  # crosspP <- apply(fdata[which(fData(es)$Entrez %in% pathways$``), ], 2, mean)

  moothaGenes <- read.table("Data/mootha.txt")
  length(moothaGenes$V1) # 447
  moothaNumbers <- which(tolower(fData(es)$gene) %in% tolower(as.character(moothaGenes$V1)))
  length(moothaNumbers) # 413
  mootha <- apply(fdata[moothaNumbers, ], 2, mean)
}

{
  p1 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(gP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Glycolysis", subtitle = "(R-MMU-70171)") +
    theme_bw(base_size = 16)
  p2 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(tcaP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("The citric acid cycle", subtitle = "(R-MMU-71403)") +
    theme_bw(base_size = 16)
  p3 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(mootha, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Oxidative phosphorylation", subtitle = "(M18264)") +
    theme_bw(base_size = 16)
  p4 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(betaP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Mitochondrial fatty acid beta-oxidation", subtitle = "(R-MMU-77289)") +
    theme_bw(base_size = 16)
  p5 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(bcaaP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("VAL, LEU and ILE degradation", subtitle = "(mmu00280)") +
    theme_bw(base_size = 16)
  p6 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(chP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Cholesterol biosynthesis", subtitle = "(R-MMU-191273)") +
    theme_bw(base_size = 16)
  p7 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(aabP, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Biosynthesis of amino acids", subtitle = "(mmu01230)") +
    theme_bw(base_size = 16)
  p8 <- pcaPlot(es.top12k, 1, 2) + aes(color=rescale(gF, to=c(-2,2))) +
    scale_color_gradientn(colours = gradientColor) +
    ggtitle("Folate metabolism", subtitle = "(mmu00790)") +
    theme_bw(base_size = 16)
}

sec_row <- plot_grid(p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     labels = c("C", "D", "E", "F"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1, 1), nrow = 1)

thd_row <- plot_grid(p5 + theme(legend.position="none"),
                     p6 + theme(legend.position="none"),
                     p7 + theme(legend.position="none"),
                     p8 + theme(legend.position="none"),
                     labels = c("G", "H", "I", "J"),
                     hjust = -2,
                     rel_widths = c(1, 1, 1, 1), nrow = 1)

q <- plot_grid(sec_row, thd_row, rel_heights = c(1, 1), ncol = 1 )





