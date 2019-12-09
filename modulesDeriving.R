# source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
source("utils.R")
show_colnames <- F
cluster_cols <- F
network <- kegg.mouse.network

solver.gmwcs <- gatom::sgmwcs.solver("sgmwcs",
                                     nthreads = detectCores(), timeLimit=600,
                                     nodes.group.by = NULL,
                                     edges.group.by = "origin",
                                     group.only.positive = T)

solver.sb <- sgmwcs.batchSolver("sgmwcs-slurm-batch",
                                nthreads=4,
                                edges.group.by = "origin", # "gene"
                                nodes.group.by = NULL,
                                group.only.positive = T,
                                c.size = 50,
                                timeLimit = 300)



# setwd(...)
# load(...)


library(org.Mm.eg.db)
entrez <- mapIds(org.Mm.eg.db,
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 keys=as.character(rownames(expTable)))
length(entrez)

rownames(expTable) <- entrez
t <- expTable
head(t)



# check spelling: entrez OR Entrez
rownames(es.top12k) <- fData(es.top12k)$entrez
head(rownames(exprs(es.top12k)))
fData(es.top12k)$entrez <- NULL
t <- exprs(es.top12k)
head(t)
dim(t)

library(GAM.networks)
data("kegg.mouse.network")
# removing HCO3-
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.x == "C00288"), ]
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.y == "C00288"), ]
network <- kegg.mouse.network



#### Preparing gene expression
# load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
gene.exprs2 <- t
gene.exprs <- t
gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$rxn2gene$gene, ]
zeroSDgenes <- (apply(gene.exprs, 1, sd, na.rm=T) == 0)
gene.exprs[zeroSDgenes,] <- gene.exprs[zeroSDgenes,] + rnorm(sum(zeroSDgenes) * ncol(gene.exprs), sd = 0.1)

# reps <- es.top12k$sample
reps <- colnames(gene.exprs)
# reps <- gsub("MOCK_.*", "MOCK", reps)

gene.cor <- cor(t(gene.exprs), use="pairwise.complete.obs")
# dim(gene.cor) # 1855 1855
gene.cor.dist <- as.dist(1 - gene.cor)
# length(gene.cor.dist)
gene.cor.dist.m <- as.matrix(gene.cor.dist)
# dim(gene.cor.dist.m) # 1855 1855
# View(head(gene.cor.dist.m)) # symmetric across diagonal


####### Initial clustering

gK <- 32
reorder <- sample(gK)
set.seed(42)
gene.pam <- pam(gene.cor.dist, k=gK)

pheatmap(
  normalize.rows(gene.exprs[gene.pam$medoids, ]),
  cluster_rows=F, cluster_cols=F,
  show_rownames=F, show_colnames=show_colnames)

##############


work.dir <- "~/Desktop/"

k <- 1
revs <- list()
curCenters <- gene.exprs[gene.pam$medoids,]
# dim(curCenters) # 16 301
# head(rownames(curCenters))

curCenters <- avearrays(curCenters, ID=reps)
curCenters <- curCenters[, reps]
base <- 0.4
dir.create(work.dir, showWarnings = F)

while (T) {
  while (T) {
    rev <- new.env()
    gK1 <- nrow(curCenters)
    rev$centers.pos <- matrix(nrow=gK1,
                              ncol=ncol(gene.exprs),
                              dimnames = list(
                                paste0("c.pos", seq_len(gK1)),
                                colnames(gene.exprs)))
    rev$centers.all <- matrix(nrow=gK1,
                              ncol=ncol(gene.exprs),
                              dimnames = list(
                                paste0("c.all", seq_len(gK1)),
                                colnames(gene.exprs)))
    rev$modules <- list()
    dist.to.centers <- 1-cor(t(curCenters), y=t(gene.exprs))
    # dim(dist.to.centers) # 16 1855
    # View(head(dist.to.centers))
    dist.to.centers[dist.to.centers < 1e-10] <- 0
    # View(dist.to.centers)
    idxs <- seq_len(gK1)
    # i = 1
    nets <- lapply(idxs, function(i) {
      messagef("Processing cluster center %s", i)
      minOther <- pmin(apply(dist.to.centers[-i, ], 2, min), base)
      score <- log2(minOther) - log2(dist.to.centers[i, ])
      # length(score) # 1855
      messagef("Number of genes scored: %s", length(score))
      # View(score)
      # removing infinities
      score[score == Inf] <- 0
      score <- pmax(score, -1000)
      es.re.scored <- makeENetworkWithScore(score, 0, network)
      es.re.scored$subnet.scored
    })
    cat("Calling solver.sb(nets)\n")
    ms <- solver.sb(nets)
    cat("Done: solver.sb(nets)\n")
    rev$modules <- ms
    for (i in idxs) {
      module <- ms[[i]]
      # capture.output(processModule(module, work.dir, sprintf("%s.gene.%s.refined.b%s", tag, i, base)))
      center.pos <- if (ulength(E(module)[score > 0]$origin) >= 3) {
        getCenter(gene.exprs, unique(E(module)[score > 0]$origin))
      } else {
        curCenters[i, ]
      }
      center.all <- if (ulength(E(module)$origin) >= 3) {
        getCenter(gene.exprs, unique(E(module)$origin))
      } else {
        curCenters[i, ]
      }
      rev$centers.pos[i, ] <- center.pos
      rev$centers.all[i, ] <- center.all
    }
    rev$centers.pos <- avearrays(rev$centers.pos, ID=reps)[, reps]
    rev$centers.all <- avearrays(rev$centers.all, ID=reps)[, reps]
    heatmapTable <- rbind(curCenters, rev$centers.pos)[rbind(
      seq_len(gK1),
      seq_len(gK1) + gK1), ]
    pheatmap(
      normalize.rows(heatmapTable),
      cluster_rows=F, cluster_cols=F,
      show_rownames=T, show_colnames=show_colnames)
    revs[[k]] <- rev
    revsToCheck <- revs[
      sapply(revs[seq_len(k-1)], function(rev) nrow(rev$centers.pos)) == nrow(rev$centers.pos)]
    diff <- max(abs(rev$centers.pos - curCenters))
    if (length(revsToCheck) > 0) {
      diff <- min(sapply(revsToCheck, function(prevRev) max(abs(rev$centers.pos - prevRev$centers.pos))))
    }
    curCenters <- rev$centers.pos
    messagef("Max diff: %s", diff)
    if (diff < 0.01) {
      break
    }
    k = k + 1
  }
  m.sizes <- sapply(revs[[k]]$modules, function(m) ulength(E(m)$origin))
  m.sizes

  if (all(m.sizes >= 5)) {
    m.sizes <- sapply(revs[[k]]$modules,
                      purrr::compose(diameter, purrr::partial(dualGraph, what="gene")))
    if (all(m.sizes >= 4)) {
      break
    }
  }

  bad <- which(m.sizes == min(m.sizes))
  centersCor <- cor(t(curCenters))
  diag(centersCor) <- NA
  toRemove <- bad[which.max(apply(centersCor, 1, max, na.rm=T)[bad])]
  curCenters <- curCenters[-toRemove, ]
  m.sizes <- m.sizes[-toRemove]
  curCenters <- curCenters[order(m.sizes, decreasing = T),]
}


# session::save.session(file="/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32_DC/session.RDa")
session::restore.session("/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/session.RDa")

# save(revs, file=sprintf("./work/gam-cluster%s.revs.rda", version))
#===================================================================================
work.dir <- "/home/octopus/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32_DC/"

# dim(kegg.mouse.network$graph.raw)
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.x == "C00288"), ]
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.y == "C00288"), ]

network <- kegg.mouse.network
# network <- kegg.human.network

curRev <- revs[[k]]
#curRev.c4 <- curRev
tag <- "m"

dir <- work.dir

for (i in seq_along(curRev$modules)) {
  message(i)
  net1 <- curRev$modules[[i]]
  V(net1)$score <- -1e-2
  m1 <- solver.gmwcs(net1)

  processModule(m1, work.dir, s=1, sprintf("%s.%s", tag, i))

  if (F && ulength(E(m1)$origin) > 100) {
    messagef("Processing cluster center %s", i)

    minOther <- pmin(apply(dist.to.centers[-i, ], 2, min), base)
    score <- log2(minOther) - log2(dist.to.centers[i, ])
    sum(score > 0)
    # removing infinities
    score[score == Inf] <- 0
    score <- pmax(score, -1000)
    es.re.scored <- makeENetworkWithScore(score, 0, network)
    m2 <- findModule(es.re.scored, solver.gmwcs)
    V(m2)$score <- -1e-3
    m2 <- solver.gmwcs(m2)

    processModule(m2, work.dir, sprintf("%s.%s.small", tag, i))

    gs <- unique(E(m2)[order(score)]$origin)
    heatmap <- gene.exprs[gs, , drop=F]
    rownames(heatmap) <- reflink[match(rownames(heatmap), Entrez), symbol]
    pheatmap(
      normalize.rows(heatmap),
      cluster_rows=F, cluster_cols=F,
      file=sprintf("%s/%s.%s.small.genes.png", work.dir, tag, i), width=12, height=8,
      show_rownames=T, show_colnames=T,
      annotation = annotation)
  }
}

# to get it go to temp.R file
orderingN <- order(ordering)
orderingN <- order(nuok)
ord <- order(colnames(expTable))
ord <- c(5,6,3,4,1,2,
         11,12,9,10,8,7,
         17,18,15,16,13,14)
annotation$organismTreatment <- ordered(annotation$organismTreatment,
                                        levels = c("none",
                                                   "IN. 2ug LPS 3 days",
                                                   "IN. 2ug LPS 6 days"))


ord <- order(colnames(curRev$centers.pos))


out <- pheatmap(
  normalize.rows(curRev$centers.pos)[, ord],
  cluster_rows=F, cluster_cols=F,
  file=sprintf("%s/%s.centers.pdf", work.dir, tag), width=25, height=10, # 12, 8
  show_rownames=T, show_colnames=T) #,
# annotation = annotation,
# annotation_colors = annotation_colors)

colnames(curRev$centers.pos[, out$tree_col[["order"]]])
col_in_order <- colnames(curRev$centers.pos[, out$tree_col[["order"]]])

#,
# legend(legend = unique(c(annotation$age,
#                        annotation$ischemic_event,
#                        annotation$`ischemic_time(days)`)))
# annotation_colors = annotation_colors)


# BUILD HEATMAPS FOR _metabolites_
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))
rma_met <- readxl::read_excel("Heart_metabolites_all_detected_571.xlsx", col_names = T)
metss <- rma_met[-1, c(1, 3, 4:35)]
colnames(metss)[c(1, 2)] <- c("name", "HMDB")

# === if we want RAW data======================
f1 <- readxl::read_excel("/home/octopus/Dropbox/GAM_RMA/H_Metabolomics_raw_data.xlsx",
                         sheet = 2, col_names = T)
f2 <- f1[6:nrow(f1), 16:ncol(f1)]
f2 <- data.matrix(f2)
# View(f2)
es <- Biobase::ExpressionSet(as.matrix(f2))
# View(es)
exprs(es) <- log2(exprs(es) + 1)
exprs(es) <- limma::normalizeBetweenArrays(exprs(es), method="quantile")
metss[, 3:ncol(metss)] <- exprs(es)
# View(metss)
# str(metss)

# === if we want QUANT NORM data ===============
f2 <- metss[, 3:ncol(metss)]
es <- ExpressionSet(as.matrix(f2))
exprs(es) <- normalizeBetweenArrays(exprs(es), method="quantile")
metss[, 3:ncol(metss)] <- exprs(es)

for (i in seq_along(curRev$modules)) {
  print(i)
  ms_KEGG <- as_ids(unique(V(curRev$modules[[i]])[order(score)]))
  print(list("ms_KEGG:", length(ms_KEGG), ms_KEGG))

  # =====================================================
  # ms_HMDB <- list()
  # for(e in seq_along(ms_KEGG)){
  #   ms_HMDB[[e]] <- met.kegg.db$mapFrom$HMDB$HMDB[
  #     which(met.kegg.db$mapFrom$HMDB$metabolite == ms_KEGG[e])
  #     ]
  # }
  # =====================================================

  # get all HMDB IDs for each KEGG ID as list
  ms_HMDB_pre <- list()
  for(e in seq_along(ms_KEGG)){
    ms_HMDB_pre[[e]] <- met.kegg.db$mapFrom$HMDB$HMDB[
      which(met.kegg.db$mapFrom$HMDB$metabolite == ms_KEGG[e])
      ]
  }

  # left only thouse in data (still list!)
  ms_HMDB_in <- lapply(ms_HMDB_pre,
                       function(x) metss[which(metss$HMDB %in% x), ])

  # ... and create a vector with its HMDB IDs
  ms_HMDB_ids_in <- lapply(ms_HMDB_pre,
                           function(x) x[which(x %in% metss$HMDB)])

  # if there is more than one row in data for one KEGG ID
  # leave one row with max sum of elements in it
  ms_HMDB_in_max <- lapply(ms_HMDB_in,
                           function(x) if(dim(x)[1] > 1) {
                             x[which.max(apply(x[, -c(1, 2)], 1, sum)), ]
                           } else {x})
  # and leave this ID in list(ms_HMDB_ids_in)
  ms_HMDB_ids_in_max <- list()
  for(j in seq_along(ms_HMDB_ids_in)) {
    if (dim(ms_HMDB_in[[j]])[1] > 1) {
      ms_HMDB_ids_in_max[[j]] <-
        ms_HMDB_ids_in[[j]][which.max(apply(ms_HMDB_in[[j]][, -c(1, 2)], 1, sum))]
    } else {ms_HMDB_ids_in_max[[j]] <- ms_HMDB_ids_in[[j]]}}

  # if one of the elemnts from list(ms_HMDB_ids_in_max) is dublicated
  # remove its element from(ms_HMDB_in_max)
  elems_to_delete <- c()
  for(k in seq_along(ms_HMDB_ids_in_max)){
    if(length(ms_HMDB_ids_in_max[[k]]) > 0){
      if(ms_HMDB_ids_in_max[[k]] %in% unlist(ms_HMDB_ids_in_max[-k]))
        elems_to_delete <- c(elems_to_delete, k)} else {next}}

  print(elems_to_delete)
  if(length(elems_to_delete > 0)){
    ms_HMDB_in_max <- ms_HMDB_in_max[-elems_to_delete]
    ms_HMDB_ids_in_max <- ms_HMDB_ids_in_max[-elems_to_delete]}

  # make normal dataframe from list for heatmap
  heatmap <- do.call(rbind, ms_HMDB_in_max)[, -c(1, 2)]
  # heatmap[, 1:8] <- rowMeans(heatmap[, 1:8])
  # heatmap[, 9:16] <- rowMeans(heatmap[, 9:16])
  # heatmap[, 17:24] <- rowMeans(heatmap[, 17:24])
  # heatmap[, 25:32] <- rowMeans(heatmap[, 25:32])
  # View(heatmap)
  rownames(heatmap) <- unname(unlist(do.call(rbind, ms_HMDB_in_max)[, 1]))

  # ====================================================
  # indxs <- which(metss$HMDB %in% unlist(ms_HMDB))
  # length(indxs)
  # heatmap <- metss[indxs, -c(1, 2)]
  # rownames(heatmap) <- unname(unlist(metss[indxs,  1]))
  # ====================================================

  pheatmap(
    normalize.rows(heatmap),
    cluster_rows=T, cluster_cols=F, # ~
    file=sprintf("%s/%s.%s.mets_max0.png", work.dir, tag, i), width=12, height=8,
    show_rownames=T, show_colnames=T)
}


ord <- order(repsss)

# BUILD HEATMAPS FOR _genes_
for (i in seq_along(curRev$modules)) {
  gs <- unique(E(curRev$modules[[i]])[order(score)]$origin)
  heatmap <- gene.exprs[gs, , drop=F]
  rownames(heatmap) <- reflink[match(rownames(heatmap), Entrez), symbol]

  # colnames(heatmap) <- gsub("WT", "WT_", colnames(heatmap))
  # colnames(heatmap) <- gsub("KO", "ERRa KO_", colnames(heatmap))
  # colnames(heatmap) <- gsub("KI", "ERBB2 KI_", colnames(heatmap))
  # colnames(heatmap) <- gsub("NE", "KI:KO_", colnames(heatmap))

  df <- normalize.rows(heatmap)
  # ord <- order(apply(df, 2, mean)) # mean

  pheatmap(
    # df,
    # df[, orderingN],
    df[, ord],
    cluster_rows=F, cluster_cols=F,
    file=sprintf("%s/%s.%s.genes.png", work.dir, tag, i), width=20, height=10, # 12, 8
    # file=sprintf("%s/MmS_%s.%s.genes_HC.png", work.dir, tag, i), width=20, height=10, # 12, 8
    show_rownames=T, show_colnames=T)
  # annotation = annotation,
  # annotation_colors = annotation_colors)
}

annotation_colors <- list()
annotation_colors$geneVar <- setNames(tol21rainbow[c(15, 16)], levels(annotation$geneVar))
annotation_colors$organ <- setNames(brewer.pal(8, "Set1"), levels(annotation$organ))
annotation_colors$cell_type <- setNames(tol21rainbow[c(11, 14, 13, 12)], levels(annotation$cell_type))




for (i in seq_along(curRev$modules)) {
  t <- get.edge.attributes(curRev$modules[[i]])[, c("origin", "symbol", "score")]
  t <- t[!duplicated(t$origin), ]
  colnames(t)[1] <- "Entrez"
  t$cor <- cor(curRev$centers.pos[i, ], t(gene.exprs[t$Entrez, ]))[1,]
  t <- t[order(t$cor, decreasing = T),]

  write.tsv(t, file=sprintf("%s/%s.%s.genes.tsv", work.dir, tag, i))

  m1 <- curRev$modules[[i]]
  net <- nets[[i]]

  notInModule <- data.table(Entrez=setdiff(E(net)$origin, E(m1)$origin))
  notInModule[, score := E(net)[match(Entrez, origin)]$score]
  notInModule[, symbol := E(net)[match(Entrez, origin)]$symbol]
  notInModule[, cor := cor(curRev$centers.pos[i, ],
                           t(gene.exprs[Entrez, ]))[1,]]
  notInModule <- notInModule[order(cor, decreasing=T), ]
  write.tsv(notInModule[score > 0], file=sprintf("%s/%s.%s.notInModule.genes.tsv", work.dir, tag, i))

}

####### --------------   org.Hs.eg.db  -----------------------------
for (i in seq_along(curRev$modules)) {
  notInModule <- data.table(Entrez=rownames(gene.exprs2))
  symbols <- mapIds(org.Mm.eg.db,
                    keys=notInModule$Entrez,
                    column="SYMBOL",
                    keytype="ENTREZID")
  notInModule[, symbol := unname(symbols)]
  notInModule[, cor := cor(curRev$centers.pos[i, ],
                           t(gene.exprs2[Entrez, ]))[1,]]
  notInModule <- notInModule[order(cor, decreasing=T), ]
  write.tsv(notInModule[1:300],
            file=sprintf("%s/%s.%s.complete.genes.tsv", work.dir, tag, i))

}
