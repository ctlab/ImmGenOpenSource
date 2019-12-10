source("utils.R")
library(reactome.db)
library(KEGGREST)
network <- kegg.mouse.network
module.dir <- "Data"
load(file="Data/243_es.top12k.Rda")
# session::restore.session("Data/session.RDa")

universe <- rownames(fData(es.top12k))
str(universe)

str(reactome.db)
reactomepath <- na.omit(AnnotationDbi::select(reactome.db, universe, "PATHID", "ENTREZID"))
reactomepath <- split(reactomepath$ENTREZID, reactomepath$PATHID)
head(reactomepath)

m_files <- list.files(module.dir, "m\\.[0-9]+\\.genes\\.tsv", full.names = T)
m_filenames <- list.files(module.dir, "m\\.[0-9]+\\.genes\\.tsv")

keggmodule <- keggLink("mmu", "module")
keggmodule <- gsub("mmu:", "", keggmodule)
names(keggmodule) <- gsub("md:", "", names(keggmodule))
keggmodule <- split(keggmodule, names(keggmodule))
# keggmodule <- lapply(listoflists, unname)
# View(keggmodule)

# keggmdnames <- KEGGREST::keggList("module", "mmu") # 404 after September, 2019
keggmdnames <- KEGGREST::keggList("module")
keggmd2name <- as.data.table(keggmdnames, keep.rownames=T)
keggmd2name$rn <- gsub("md:", "", keggmd2name$rn)
setnames(keggmd2name, c("rn","keggmdnames"), c("PATHID","PATHNAME"))
keggmd2name$PATHID <- paste0("mmu_", keggmd2name$PATHID)
# View(keggmd2name)

keggpathway <- keggLink("mmu", "pathway")
keggpathway <- gsub("mmu:", "", keggpathway)
names(keggpathway) <- gsub("path:", "", names(keggpathway))
keggpathway <- split(keggpathway, names(keggpathway))
keggpathway <- lapply(keggpathway, unname)
keggpathnames <- keggList("pathway", "mmu")

keggpath2name <- as.data.table(keggpathnames, keep.rownames=T)
keggpath2name$rn <- gsub("path:", "", keggpath2name$rn)
keggpath2name$keggpathnames <- gsub(" - Mus musculus \\(mouse\\)", "",
                                    keggpath2name$keggpathnames)
setnames(keggpath2name, c("rn","keggpathnames"), c("PATHID","PATHNAME"))

reactomepathway2name <- as.data.table(na.omit(
  AnnotationDbi::select(reactome.db,
                        names(reactomepath),
                        c("PATHNAME"), 'PATHID')))
# Remove organism prefix
reactomepathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]

# combine kegg modules and pathways with reactome data
pathways <- c(reactomepath, keggmodule, keggpathway)
pathways <- pathways[sapply(pathways, length) >= 10]
#
pathways$`5991024` <- NULL
pathways$`R-MMU-1430728` <- NULL # Metabolism
pathways$`mmu01100` <- NULL # Metabolic pathways
pathways$`mmu01200` <- NULL # Carbon metabolism
# pathways$`mmu01230` <- NULL # Biosynthesis of amino acids
# pathways$`R-MMU-71387` <- NULL # Metabolism of carbohydrates

pathway2name <- do.call("rbind", list(reactomepathway2name,
                                      keggmd2name,
                                      keggpath2name))

gseaReactome <- function(genes) {
  overlaps <- data.frame(
    q=sapply(sapply(pathways, intersect, genes), length),
    m=sapply(pathways, length),
    n=length(universe)-sapply(pathways, length),
    k=length(genes))

  igenes <- sapply(pathways, intersect, genes)

  for(i in seq_along(igenes)){
    overlaps[i, 5] <- paste(igenes[[i]], collapse=" ")}

  # q-1 because we want probability of having >=q white balls
  pathways.pvals <- with(overlaps,
                         mapply(phyper, q-1, m, n, k, lower.tail = FALSE))
  res <- data.table(PATHID=names(pathways),
                    pval=pathways.pvals,
                    k=overlaps$q,
                    K=overlaps$m,
                    genes = overlaps$V5)

  res[, padj := p.adjust(pval, method="BH")]

  res <- merge(res, pathway2name, by="PATHID")
  res <- res[order(pval),]
  res <- res[, c(1, 2, 3, 4, 6, 7, 5)]
}

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
# dplyr::glimpse(org.Mm.eg.gatom.anno)

for (i in 1:length(m_files)){
  file <- data.table::fread(m_files[i])
  name <- basename(m_files[i])
  name <- sapply(strsplit(name, ".", fixed = T),"[[", 2)
  out <- gseaReactome(file$Entrez)

  out <- out[padj < 0.05,]

  pz <- sapply(out$genes, function(x) strsplit(x, " "))
  new_pz <- vector("list", length = nrow(out))
  for(e in seq_along(pz) ){
    for(j in seq_along(pz[[e]]) ){
      new_pz[[e]] <- append(new_pz[[e]],
                            org.Mm.eg.gatom.anno$genes$symbol[which(org.Mm.eg.gatom.anno$genes$gene
                                                                    == pz[[e]] [[j]] )] ) } }

  for(a in seq_along(new_pz)){
    out[a, 7] <- paste(new_pz[[a]], collapse=" ")}

  write.tsv(out,
            file=sprintf("%s/m.%s.pathways_mod.tsv", module.dir, name))
  print(basename(m_files[i]))
}



