source("utils.R")
load("Data/243_es.top12k.Rda")
work.dir <- "Data"
solver.gmwcs <- gatom::sgmwcs.solver("sgmwcs",
                                     nthreads=detectCores(), 
                                     timeLimit=600,
                                     nodes.group.by=NULL,
                                     edges.group.by="origin",
                                     group.only.positive = T)
solver.sb <- sgmwcs.batchSolver("sgmwcs-slurm-batch",
                                nthreads=4,
                                edges.group.by="origin", # "gene"
                                nodes.group.by=NULL,
                                group.only.positive=T,
                                c.size=50,
                                timeLimit=300)

rownames(es.top12k) <- fData(es.top12k)$entrez
head(rownames(exprs(es.top12k)))
fData(es.top12k)$entrez <- NULL
t <- exprs(es.top12k)
head(t)
dim(t)

kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.x == "C00288"), ]
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.y == "C00288"), ]
network <- kegg.mouse.network
# network <- kegg.human.network

### Preparing gene expression
gene.exprs2 <- t
gene.exprs <- t
gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$rxn2gene$gene, ]
zeroSDgenes <- (apply(gene.exprs, 1, sd, na.rm=T) == 0)
gene.exprs[zeroSDgenes,] <- gene.exprs[zeroSDgenes,] + rnorm(sum(zeroSDgenes) * ncol(gene.exprs), sd = 0.1)

reps <- colnames(gene.exprs)

gene.cor <- cor(t(gene.exprs), use="pairwise.complete.obs")
gene.cor.dist <- as.dist(1 - gene.cor)
gene.cor.dist.m <- as.matrix(gene.cor.dist)

### Initial clustering
gK <- 32
reorder <- sample(gK)
set.seed(42)
gene.pam <- pam(gene.cor.dist, k=gK)

pheatmap(
  normalize.rows(gene.exprs[gene.pam$medoids, ]),
  cluster_rows=F, cluster_cols=F,
  show_rownames=F, show_colnames=F)

k <- 1
revs <- list()
curCenters <- gene.exprs[gene.pam$medoids,]

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
    dist.to.centers[dist.to.centers < 1e-10] <- 0
    idxs <- seq_len(gK1)
    nets <- lapply(idxs, function(i) {
      messagef("Processing cluster center %s", i)
      minOther <- pmin(apply(dist.to.centers[-i, ], 2, min), base)
      score <- log2(minOther) - log2(dist.to.centers[i, ])
      messagef("Number of genes scored: %s", length(score))
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
      show_rownames=T, show_colnames=F)
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

# session::save.session(file=sprintf("%s/session.Rda", work.dir))
# session::restore.session(sprintf("%s/session.Rda", work.dir))

curRev <- revs[[k]]
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

out <- pheatmap(
  normalize.rows(curRev$centers.pos),
  cluster_rows=F, cluster_cols=F,
  file=sprintf("%s/%s.centers.pdf", work.dir, tag), width=25, height=10,
  show_rownames=T, show_colnames=T) 

### heatmap for genes
for (i in seq_along(curRev$modules)) {
  gs <- unique(E(curRev$modules[[i]])[order(score)]$origin)
  heatmap <- gene.exprs[gs, , drop=F]
  rownames(heatmap) <- reflink[match(rownames(heatmap), Entrez), symbol]

  df <- normalize.rows(heatmap)

  pheatmap(
    df,
    cluster_rows=F, cluster_cols=F,
    file=sprintf("%s/%s.%s.genes.png", work.dir, tag, i), width=20, height=10,
    show_rownames=T, show_colnames=T)
}

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

### use org.Hs.eg.db in case of human samples
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
