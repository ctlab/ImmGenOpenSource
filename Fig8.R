library(pheatmap)
library(RColorBrewer)

load("~/Documents/immGen/final_objects/337_es.total.Rda")

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



# genes
fData(es)$gene[grep("Nos|nos", fData(es)$gene)]
gnsh <- c("Il1b", "Il6", "Acod1", "Arg1", "Tnf", "Cxcl2", "Cxcl10", "Ptgs2",
          "Retnla", "Ccl22", "H2-Ab1", "H2-Aa", "Mrc1", "Clec10a", "Nos2")

ttt <- exprs(es)[which(fData(es)$gene %in% gnsh), which(es$treatment %in%
                                                          c("Control (C. albicans)",
                                                            "8h post-infection (C. albicans)",
                                                            "48h post-infection (C. albicans)"))]
dim(ttt)
colnames(ttt)
ttt <- ttt[, calbList]

cr <- colorRampPalette(c("#0000ff","#ffffff","#ff0000"))

out <- pheatmap(
  rUtils::normalize.rows(ttt),
  cluster_rows=T, cluster_cols=F,
  file="~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig8/genesHeatmap.pdf", width=5, height=5, # 12, 8
  show_rownames=T, show_colnames=T,
  color = cr(1000), border_color = "grey60") #,




# modules
# load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/session.RDa")
session::restore.session("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/session.RDa")

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
  file="~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig8/modulesHeatmap.pdf", width=5, height=4.25,
  show_rownames=T, show_colnames=T,
  color = cr(1000))
  # color = rev(brewer.pal(9,"RdYlBu")))




# pathways
load("~/Documents/immGen/pathways/pathway2name.rda")
load("~/Documents/immGen/pathways/pathways.rda")

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
  file="~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig8/pthwsHeatmap.pdf", width=6, height=4, # 12, 8
  show_rownames=T, show_colnames=T,
  color = cr(1000), border_color = "grey60")
  # color = rev(brewer.pal(9,"RdYlBu")))



# DE
### Loading count data
countData <- read.table("~/Documents/immGen/OSMNP_unnormalized_genes_count_10_3_18.count_table", sep = "\t", header = T)
rownames(countData) <- countData$gene_symbol
countData <- countData[,-1]
dim(countData) # 52997   414

### Loading meta data
load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/es.top12k.Rda")
# load("~/Documents/immGen/final_objects/337_es.top12k.Rda") # the same but not almost correct metaSample anno
anno <- Biobase::pData(es.top12k)
unique(anno$metaSample)
View(anno)
# cellTypes
anno$cellTypes <- anno$metaSample
anno$cellTypes[anno$metaSample %in% c("MF_tissue_resident",
                                      "MF_embr_E6.0E8.0_and_Alv_and_SPM",
                                      "MG",
                                      "MF_adipose_tissue")] <- "MF"
anno$cellTypes[anno$metaSample %in% c("DC_tissue_resident",
                                      "pDC_BM_Sp",
                                      "DC_mig")] <- "DC"
unique(anno$cellTypes)

### Filter and reorder count table
countData <- countData[, rownames(anno)]
dim(countData)
head(colnames(countData))
head(rownames(anno))

### Filter rarely expressed genes
y <- edgeR::DGEList(counts = countData)
keep <- rowSums(edgeR::cpm(y) > 1) >= 3
y <- y[keep,]
dim(y)
y$samples$lib.size = colSums(y$counts) # reset lib sizes

### Normalisation and log
ytmm <- edgeR::calcNormFactors(y, method = "TMM")
# limma::plotMDS(ytmm, col=anno$Batch, pch=16) # metric - log2 FC-s between samples

### Gene filtering
entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                keytype = "SYMBOL",
                                column = "ENTREZID",
                                keys=as.character(rownames(ytmm)))
unique(entrez[which(duplicated(entrez))]) # NA
length(which(is.na(entrez))) #

ytmm <- ytmm[-which(is.na(entrez)), ]
ytmm <- ytmm[head(order(apply(ytmm, 1, mean), decreasing = T), 12000), ]
dim(ytmm)
head(rownames(ytmm))

### log and de
design <- model.matrix(~0+cellTypes, data = anno); design
design <- model.matrix(~0+metaSample, data = anno); design
# colnames(design) <- levels(TS)
cs <- unique(colnames(design)); cs
dplyr::all_equal(cs, colnames(design))

library(limma)
v <- voom(ytmm, design, plot = TRUE)
fit <- lmFit(v, design)

des <- list()
cntrsts <- c()

for(i in seq_along(cs)){
  for(j in seq_along(cs)){
    if(i==j)
      next
    cntrsts <- c(cntrsts, paste(cs[[i]], "-", cs[[j]]))
  }
}

for(i in seq_along(cntrsts)){
  crrnt <- cntrsts[i]
  fit2 <- contrasts.fit(fit, makeContrasts(vs=crrnt,
                                           levels=design))
  fit2 <- eBayes(fit2)
  de <- topTable(fit2, adjust.method="BH", number=Inf)
  tag <- sprintf("%s.vs.%s",
                 strsplit(crrnt, " - ")[[1]][2],
                 strsplit(crrnt, " - ")[[1]][1])
  des[[tag]] <- de

  rUtils::write.tsv(de,
                    file=
                      sprintf("~/Documents/immGen/DE/20191105_metaSample/cellTypes/%s.tsv",
                              tag))
}






# gatom / GAM
path <- "~/Documents/immGen/DE/20191105_metaSample/metaSamples/need/"
path <- "~/Documents/immGen/DE/20191105_metaSample/cellTypes/"
path <- "~/Documents/immGen/DE/20191105_metaSample/metaSamples/Mo_vs_allMFs/"
filenames <- list.files(path, pattern="*.tsv", full.names=TRUE)
list_of_de2 <- lapply(filenames, readr::read_tsv)

k.gene <- c(85) #,100

nms <- gsub(".tsv", "",
            gsub("metaSample", "",
                 list.files(path, pattern="*.tsv")))
nms

library(dplyr)
library(gatom)
library(data.table)
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
org.Mm.eg.gatom.anno$mapFrom$Symbol <- rbind(org.Mm.eg.gatom.anno$mapFrom$Symbol, data.table(gene="16365", Symbol="Irg1"))
data.table::setkeyv(org.Mm.eg.gatom.anno$mapFrom$Symbol, "Symbol")

# i = 3
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















