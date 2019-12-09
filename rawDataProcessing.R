### Loading count data
countData <- read.table("~/Documents/immGen/OSMNP_unnormalized_genes_count_10_3_18.count_table", sep = "\t", header = T)
dim(countData)
View(head(countData))
which(duplicated(countData$gene_symbol)==T) # no duplicated rows
rownames(countData) <- countData$gene_symbol
countData <- countData[,-1]
dim(countData) # 52997   414



### Loading meta data
anno <- read.csv("~/Documents/immGen/Official_OSMNP_metadata_QC_passes6.19.18.csv",
                 header = TRUE, stringsAsFactors = FALSE)
dim(anno) # 417  20
View(anno)
# fixing typos, etc
anno$Organ[anno$Organ %in% c("Peritoenal Cavity", "Peritoneal")] <- "Peritoneal Cavity"
anno$Organ[which(anno$Organ == "yolk sac")] <- "Yolk sac"
anno$Organ[which(anno$Organ == "Epithelial Cell")] <- "Epithelium"
anno$Organ[which(anno$Organ == "Lymph Node (Mesenteric)")] <- "Mesenteric LN"
anno$Organ[which(anno$Organ == "Mesenteric Fat")] <- "Mesenteric fat"
rownames(anno) <- anno$ImmGen.Nomenclature
View(head(anno))



### Filtering samples in anno table
anno <- anno[-which(!rownames(anno) %in% colnames(countData)), ]
anno <- anno[-which(rownames(anno) %in% c("MF.FC.1", "MF.FC.3")), ] # contamination with fibroblasts
View(anno)
PC_to_keep <- c("MF.226pIIp.PC.1",
                "MF.226pIIp.PC.3",
                "MF.226pIIp.PC.4",
                "MF.226pIIp480lo.PC.1",
                "MF.226pIIp480lo.PC.2", # SPMs from Ki Wook, Randolph Lab
                "MF.PC.01", # Jigar Desai from Lionakis lab, ICAM2+F4/80+
                "MF.PC.03", # Ivan Dzhagalov from Dzhagalov lab
                "MF.PC.07", # Geetika Bajpai from Lavine lab
                "MF.PC.09") # Karine Crozat from Dalod lab
anno <- anno[-setdiff(
  which(anno$Organ == "Peritoneal Cavity"),
  which(anno$ImmGen.Nomenclature %in% PC_to_keep)), ] # dowsampling MF.PC
View(anno)

tissue1 <- dplyr::as_tibble(readr::read_csv(
  "/home/octopus/Documents/immGen/Official_OSMNP_master_table6.19.18.csv"))
View(tissue1)

# (var1) Keeping all samples
anno$treatment <- tissue1$Organism.Treatment[match(anno$SampleName, tissue1$Sample.me)]
# tissue1$Cell.Treatment
# tissue1$Genetic.Variation

# (var2) Removing treated samples
anno <- anno[-grep("pIC.alv", anno$ImmGen.Nomenclature), ]
treated <- tissue1$Sample.me[which(tissue1$Organism.Treatment != "none")]
anno <- anno[-which(anno$SampleName %in% treated), ]

# (var3) Removing treated and embryo samples
anno <- anno[-which(anno$Organ %in% c("Yolk sac", "Cell (Embryonic Body)")), ]
anno <- anno[-grep("E[1468]{2}|Neo|pIC.alv", anno$ImmGen.Nomenclature), ]
treated <- tissue1$Sample.me[which(tissue1$Organism.Treatment != "none")]
anno <- anno[-which(anno$SampleName %in% treated), ]

# Test filtering of MG (not used)
# MG_to_keep <- c("MF.45lo.CNS.1",
#                 "MF.microglia.cerebel.CNS.1",
#                 "MF.microglia.cerebr.CNS.1",
#                 "MG.cortex.1",
#                 "MG.hippo.1",
#                 "MF.microglia.CNS.1")
# anno <- anno[-setdiff(
#   which(anno$Organ == "Brain"),
#   which(anno$ImmGen.Nomenclature %in% MG_to_keep)), ]


dim(anno) # 204 20 (for General), 243 20 (for GAM-clustering), # 337 20 (for Metab)
View(anno)
anno$treatment[anno$ImmGen.Nomenclature == "MF.pIC.alv.siglecFp.Lu.2"] <- "pIC"
anno$Cell.Family[anno$Organ == "Brain"] <- "Microglia"
anno$Sorting.Markers <- gsub("\t", "", anno$Sorting.Markers)
anno$Sorting.Markers <- gsub("\\s*", "", anno$Sorting.Markers)
anno$treatment[is.na(anno$treatment)] <- "none"
rUtils::write.tsv(anno, file = "~/Documents/immGen/final_objects/337_anno.tsv")



### Filter and reorder count table
countData <- countData[, rownames(anno)]
dim(countData)



### Filter rarely expressed genes
eObj <- edgeR::DGEList(counts = countData)
keep <- rowSums(edgeR::cpm(eObj) > 1) >= 3
eObj <- eObj[keep,]
dim(eObj)
eObj$samples$lib.size = colSums(eObj$counts) # reset lib sizes



### Normalisation and log
ytmm <- edgeR::calcNormFactors(eObj, method = "TMM")
limma::plotMDS(ytmm, col=anno$Batch, pch=16) # metric - log2 FC-s between samples
# .
# TS <- paste(anno$Cell.Family, anno$Organ, anno$Lab, anno$Batch, sep=".")
# TS <- paste(anno$Cell.Family, anno$Organ, anno$Lab, sep=".")
# TS <- factor(TS, levels=unique(TS))
# design <- model.matrix(~0+TS); design
# colnames(design) <- levels(TS)
# v <- limma::voom(ytmm, design, plot = TRUE)
# .
v <- limma::voom(ytmm, plot = TRUE)
expTable <- v$E
dim(expTable)


library(Biobase)
es <- ExpressionSet(assayData = expTable,
                    phenoData = new("AnnotatedDataFrame",
                                     data = anno))
fData(es)$gene <- rownames(expTable)


### Filtering genes
fData(es)$entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                          keytype = "SYMBOL",
                                          column = "ENTREZID",
                                          keys=as.character(fData(es)$gene))
unique(fData(es)$entrez[which(duplicated(fData(es)$entrez))]) # NA
length(which(is.na(fData(es)$entrez))) #

#
save(es, file = "~/Documents/immGen/final_objects/337_es.total.Rda")

es <- es[-which(is.na(fData(es)$entrez)), ]
es.top12k <- es[head(order(apply(exprs(es), 1, mean), decreasing = T), 12000), ]


#
save(es.top12k, file = "~/Documents/immGen/final_objects/337_es.top12k.Rda")
# save(es.top12k, file = "~/Documents/immGen/final_objects/243_es.top12k.Rda")
# save(es.top12k, file = "~/Documents/immGen/final_objects/150MF_es.top12k.Rda")
# save(es.top12k, file = "~/Documents/immGen/final_objects/91DC_es.top12k.Rda")








