path <- "~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/"
m1 <- rUtils::read.tsv(paste0(path, "m.1.genes.tsv"))
m2 <- rUtils::read.tsv(paste0(path, "m.2.genes.tsv"))
m8 <- rUtils::read.tsv(paste0(path, "m.8.genes.tsv"))
m9 <- rUtils::read.tsv(paste0(path, "m.9.genes.tsv"))

intersect(m1$symbol, m2$symbol) # 6. ~half of them will be removed

intersect(m1$symbol, m9$symbol)
intersect(m2$symbol, m9$symbol)

intersect(m1$symbol, m8$symbol)
intersect(m2$symbol, m8$symbol)

intersect(m9$symbol, m8$symbol)


#--------------------------------------------------

session::restore.session("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/session.RDa")

work.dir <- "~/Desktop/"
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.x == "C00288"), ]
kegg.mouse.network$graph.raw <- kegg.mouse.network$graph.raw[-which(
  kegg.mouse.network$graph.raw$met.y == "C00288"), ]
network <- kegg.mouse.network
curRev <- revs[[k]]
tag <- "m"
dir <- work.dir

#--------------------------------------------------

net1 <- curRev$modules[[1]]
net2 <- curRev$modules[[2]]
# as_data_frame(net1)$symbol[duplicated(as_data_frame(net1)$symbol)]

relations <- rbind(
  cbind(as_data_frame(net1), net="net1"),
  cbind(as_data_frame(net2), net="net2"))

View(relations)

relations$uni <- paste(relations$from, relations$to, relations$symbol, sep="_")
relations$uni[duplicated(relations$uni)]
View(relations)
dublRows <- which(relations$uni %in% relations$uni[duplicated(relations$uni)])
View(relations[dublRows, -which(colnames(relations)%in%c("pathway",
                                                         "from",
                                                         "to",
                                                         "symbol",
                                                         "label",
                                                         "origin"))])
remRows <- rownames(relations[dublRows, ][which(relations[dublRows, ]$score < 0), ]) # but it goesnt matter because we will
# rewrite scores, what is important is that we are removing duplicated rows
relations12pre <- relations[-as.numeric(remRows), -which(colnames(relations)=="uni")]









# average centers and recalculate correlation
curCenter12 <- apply(curCenters[c(1,2), ], 2, mean)
curCenters[1,] <- curCenter12
curCenters <- curCenters[-2,]
View(curCenters)

dist.to.centers <- 1-cor(t(curCenters), y=t(gene.exprs))
dim(dist.to.centers) # 8 1835
dist.to.centers[dist.to.centers < 1e-10] <- 0
i=1
minOther <- pmin(apply(dist.to.centers[-i, ], 2, min), base)
score <- log2(minOther) - log2(dist.to.centers[i, ])
score[score == Inf] <- 0
score <- pmax(score, -1000)


relations12pre$newScore <- unname(score[match(relations12pre$gene, names(score))])

View(relations12pre[, -which(colnames(relations)%in%c("pathway",
                                                 "from",
                                                 "to",
                                                 "symbol",
                                                 "label",
                                                 "origin"))])
View(relations12pre)

relations12 <- relations12pre
relations12$score <- relations12pre$newScore
relations12$log2FC <- relations12pre$newScore
relations12$pval <- 2^(-relations12$score)
relations12$logPval <- log(relations12$pval)
View(relations12)
relations12 <- relations12[, -which(colnames(relations12pre)%in%c(
                                                           "pval",
                                                           "log2FC",
                                                           "logPval",
                                                           "uni",
                                                           "newScore",
                                                           "net"
                                                           ))]
View(relations12)
dim(relations12) # 87, 11
duplicated(relations12)
# relations12 <- relations12[!duplicated(relations12), ]





# es.re.scored <- makeENetworkWithScore(relations12$score, 0, network)
# es.re.scored$subnet.scored











actors <- rbind(
  as_data_frame(net1, what = "vertices"),
  as_data_frame(net2, what = "vertices"))
View(actors)
duplicated(actors)
actors12 <- actors[!duplicated(actors), ]

net12 <- graph_from_data_frame(relations12, directed=F, vertices=actors12)

m12 <- solver.gmwcs(net12)

i=1.2
processModule(m12, work.dir, s=1, sprintf("%s.%s", tag, i))







View(rbind(
  curCenters,
  curRev$centers.pos))

out <- pheatmap(
  normalize.rows(curCenters),
  cluster_rows=F, cluster_cols=F,
  file=sprintf("%s/%s.centers.pdf", work.dir, tag), width=25, height=10, # 12, 8
  show_rownames=T, show_colnames=T) #,


# ------------------------------------------


i=1

# ------------------------------------------


t <- get.edge.attributes(m12)[, c("origin", "symbol", "score")]
t <- t[!duplicated(t$origin), ]
colnames(t)[1] <- "Entrez"
t$cor <- cor(curCenters[i, ], t(gene.exprs[t$Entrez, ]))[1,]
t <- t[order(t$cor, decreasing = T),]

write.tsv(t, file=sprintf("%s/%s.%s.genes.tsv", work.dir, tag, i))

m1 <- m12 # NB m12 IS NOT THE SAME AS m1 <- curRev$modules[[i]], actually net12 is the same
net <- nets[[i]]

notInModule <- data.table(Entrez=setdiff(E(net)$origin, E(m1)$origin))
notInModule[, score := E(net)[match(Entrez, origin)]$score]
notInModule[, symbol := E(net)[match(Entrez, origin)]$symbol]
notInModule[, cor := cor(curCenters[i, ],
                         t(gene.exprs[Entrez, ]))[1,]]
notInModule <- notInModule[order(cor, decreasing=T), ]
write.tsv(notInModule[score > 0], file=sprintf("%s/%s.%s.notInModule.genes.tsv", work.dir, tag, i))

# ------------------------------------------


notInModule <- data.table(Entrez=rownames(gene.exprs2))
symbols <- mapIds(org.Mm.eg.db,
                  keys=notInModule$Entrez,
                  column="SYMBOL",
                  keytype="ENTREZID")
notInModule[, symbol := unname(symbols)]
notInModule[, cor := cor(curCenters[i, ],
                         t(gene.exprs2[Entrez, ]))[1,]]
notInModule <- notInModule[order(cor, decreasing=T), ]
write.tsv(notInModule[1:300],
          file=sprintf("%s/%s.%s.complete.genes.tsv", work.dir, tag, i))

# ------------------------------------------
