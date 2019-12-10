source("utils.R")
load("Data/pathway2name.rda")
load("Data/pathways.rda")
load("Data/337_es.total.Rda")

Control_Calb <- c(
  "MF.45lo.CNS.1",
  "MF.45lo.CNS.2",
  "MF.6Gn480hi.Kd.1",
  "MF.6Gn480hi.Kd.2",
  "MF.11blop64p169p.Lv.1",
  "MF.11blop64p169p.Lv.2")
es$treatment[which(es$ImmGen.Nomenclature %in% Control_Calb)] <- "Control (C. albicans)"

calbList <- c("MF.11blop64p169p.Lv.1",
              "MF.11blop64p169p.Lv.2",
              "MF.11blop64p169p.Calb.8h.Lv.1",
              "MF.11blop64p169p.Calb.8h.Lv.2",
              "MF.11blop64p169p.Calb.48h.Lv.1",
              "MF.11blop64p169p.Calb.48h.Lv.2",
              "MF.6Gn480hi.Kd.1",
              "MF.6Gn480hi.Kd.2",
              "MF.6Gn480hi.Calb.8h.Kd.1",
              "MF.6Gn480hi.Calb.8h.Kd.2",
              "MF.6Gn480hi.Calb.48h.Kd.1",
              "MF.6Gn480hi.Calb.48h.Kd.2",
              "MF.45lo.CNS.1",
              "MF.45lo.CNS.2",
              "MF.45lo.Calb.8h.CNS.1",
              "MF.45lo.Calb.8h.CNS.2",
              "MF.45lo.Calb.48h.CNS.1",
              "MF.45lo.Calb.48h.CNS.2")

### genes
gnsh <- c("Il1b", "Il6", "Acod1", "Arg1", "Tnf", "Cxcl2", "Cxcl10", "Ptgs2",
          "Retnla", "Ccl22", "H2-Ab1", "H2-Aa", "Mrc1", "Clec10a", "Nos2")

ttt <- exprs(es)[which(fData(es)$gene %in% gnsh), calbList]

cr <- colorRampPalette(c("#0000ff","#ffffff","#ff0000"))

out <- pheatmap(
  rUtils::normalize.rows(ttt),
  cluster_rows=T, cluster_cols=F,
  file="genesHeatmap.pdf", width=5, height=5,
  show_rownames=T, show_colnames=T,
  color = cr(1000), border_color = "grey60")

### modules
session::restore.session("Data/session.RDa")

curRev <- revs[[k]]
range(gene.exprs)
expTable <- exprs(es)[, calbList]
dim(expTable)

expTable <- t(apply(expTable, 1, function (x) mosaic::zscore(x)))
eyegene <- as.data.frame(matrix(nrow=length(curRev$modules), ncol=18))

for (i in seq_along(curRev$modules)) {
  print(paste("module #", i))
  gs <- unique(E(curRev$modules[[i]])[order(score)]$origin)
  print(length(gs))
  heatmap <- expTable[which(fData(es)$entrez %in% gs), , drop=F]
  eyegene[i, ] <- apply(heatmap, 2, mean)
}
colnames(eyegene) <- colnames(expTable)
rownames(eyegene) <- paste("Module", 1:9)

out <- pheatmap(
  rUtils::normalize.rows(eyegene),
  cluster_rows=F, cluster_cols=F,
  file="modulesHeatmap.pdf", width=5, height=4.25,
  show_rownames=T, show_colnames=T,
  color = cr(1000))
  # color = rev(brewer.pal(9,"RdYlBu")))

### pathways
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
pathway2name$PATHNAME[grep("folate", pathway2name$PATHNAME)]
fP <- pathways$`mmu00790`
pathway2name$PATHID[which(pathway2name$PATHNAME == "Metabolism of folate and pterines")]
fpP <- pathways$`R-MMU-196757`

pathway2name$PATHID[which(pathway2name$PATHNAME == "Signaling by Interleukins")]
ilP <- pathways$`R-MMU-449147`
pathway2name$PATHID[which(pathway2name$PATHNAME == "Antigen processing-Cross presentation")]
mhcP <- pathways$`R-MMU-1236975`
pathway2name$PATHID[which(pathway2name$PATHNAME == "NF-kappa B signaling pathway")]
nfkbP <- pathways$`mmu04064`

PATs <- rep(NA, 18)
pathsh <- list(gP, pppP, tcaP, betaP, bcaaP, chP, aabP, fpP, ilP, mhcP, nfkbP)
pnsh <- c("Glycolysis (R-MMU-70171)", "PPP (mmu00030)", "TCA (R-MMU-71403)",
          "FA beta-oxidation (R-MMU-77289)", "BCAA degradation (mmu00280)", "Cholesterol metabolism (R-MMU-191273)",
          "Biosynthesis of amino acids (mmu01230)", "Folate metabolism (R-MMU-196757)",
          "IL-signaling (R-MMU-449147)", "Cross presentation (R-MMU-1236975)", "NF-kappaB signaling (mmu04064)")

for (i in seq_along(pathsh)) {
  print(pnsh[i])
  print(length(pathsh[[i]]))
  print(length(which(fData(es)$entrez %in% pathsh[[i]])))
  heatmap0 <- expTable[which(fData(es)$entrez %in% pathsh[[i]]), ]
  heatmap1 <- apply(heatmap0, 2, mean)
  PATs <- rbind(PATs, heatmap1)
}
# View(PATs)
PATs <- PATs[-1, ]
rownames(PATs) <- pnsh

out <- pheatmap(
  rUtils::normalize.rows(PATs),
  cluster_rows=F, cluster_cols=F,
  file="pthwsHeatmap.pdf", width=6, height=4, # 12, 8
  show_rownames=T, show_colnames=T,
  color = cr(1000), border_color = "grey60")
  # color = rev(brewer.pal(9,"RdYlBu")))

### DE
countData <- read.table("Data/OSMNP_unnormalized_genes_count_10_3_18.count_table", sep = "\t", header = T)
# <making ExpressionSet with unnormalized counts and taking top 12k genes>
es_dds <- es.top12k[, calbList]

pData(es_dds)$C.alb <- c(
 "LvContr", "LvContr", "Lv8h", "Lv8h",  "Lv48h", "Lv48h",
 "KdContr", "KdContr", "Kd8h", "Kd8h", "Kd48h", "Kd48h",  
 "CNSContr", "CNSContr", "CNS8h", "CNS8h", "CNS48h", "CNS48h")

dds <- rUtils::DESeqDataSetFromExpressionSet(es_dds, design = ~ C.alb)
dds <- dds[rowSums(counts(dds)) > 20, ]
dds <- DESeq(dds)

cs <- unique(es_dds$LPS)
des <- list()
for (i in seq_along(cs)) {
  for (j in seq_along(cs)) {
    if (i == j)
      next
    
    res <- rUtils::normalizeGeneDE(results(dds, contrast = c("LPS", cs[j], cs[i]), cooksCutoff = F))
    res <- res[head(order(baseMean, decreasing=T), 12000), ]
    res <- res[order(stat), ]
    tag <- sprintf("%s.vs.%s", cs[i], cs[j])
    des[[tag]] <- res
    
    rUtils::write.tsv(res, file=sprintf(paste0(current_path, "DESeq_%s.tsv"), tag))
  }
}

# gatom / GAM
path <- "Data"
filenames <- list.files(path, pattern="DESeq_*.tsv", full.names=TRUE)
list_of_de2 <- lapply(filenames, readr::read_tsv)

k.gene <- c(60, 70, 85, 95)

nms <- gsub(".tsv", "", gsub("metaSample", "", list.files(path, pattern="*.tsv")))

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
org.Mm.eg.gatom.anno$mapFrom$Symbol <- rbind(org.Mm.eg.gatom.anno$mapFrom$Symbol, data.table(gene="16365", Symbol="Irg1"))
data.table::setkeyv(org.Mm.eg.gatom.anno$mapFrom$Symbol, "Symbol")

ed2 <- function(list_of_de2, k.gene, nms) {

  for (i in seq_along(list_of_de2)) {

    de <- list_of_de2[[i]]
    name <- nms[i]

    print(de[de$AveExpr == 0, ])
    print(de[de$P.Value == 0, ])
    print(de[de$baseMean == 0, ])
    print(de[de$pval == 0, ])

    # de_data2 <- de
    de_data2 <- de %>%
      dplyr::select(X1, AveExpr, logFC, P.Value) %>%
      dplyr::rename(ID = X1,
                    baseMean = AveExpr,
                    log2FC = logFC,
                    pval = P.Value)

    # de_data2 <- de %>%
    #   dplyr::select(ID, baseMean, log2FC, pval)

    head(de_data2)

    g <- makeAtomGraph(network = network,
                       org.gatom.anno = org.Mm.eg.gatom.anno,
                       gene.de = de_data2,
                       met.db = met.kegg.db,
                       met.de = NULL)

    for (k in k.gene) {
      # k=85
      gs <- scoreGraph(g,
                       k.gene = k,
                       k.met = NULL)
      solve <- sgmwcs.solver("sgmwcs",
                             nthreads = 4,
                             timeLimit = 60,
                             minimize.size = TRUE)

      # solve <- sgmwcs.solver("sgmwcs", nthreads = 4, timeLimit = 30)
      # library(igraph)
      # V(gs)$score <- -0.001
      # solve <- sgmwcs.solver("sgmwcs", nthreads = 4, timeLimit = 30, minimize.size = F)

      # print("solving graph %s with k=%s..." %f% c(name, k))
      m <- solve(gs); m

      # print("printing graph %s with k=%s..." %f% c(name, k))
      # for(sss in seq_along(1:5)+75){
      sss <- 22
      m1 <- addHighlyExpressedEdges(m, gs, top=12000); m1
      m1 <- collapseAtomsIntoMetabolites(m1); m1
      saveModuleToPdf(m, n_iter = 120, force = 1e-5, seed = sss,
                      file = paste0(path, name, "_k", k, "_seed", sss, ".pdf"),
                      name = name)
      # }
      gatom::saveModuleToXgmml(module = m,
                               file = paste0(path, name, "_k", k, ".xgmml"))
    }
  }
}

ed2(list_of_de2, k.gene, nms)