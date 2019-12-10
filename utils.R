Packages <- c("BiocManager", "BioNet", "Biobase", 
              "devtools", "parallel", "logging", 
              "cluster", "data.table", "plyr", 
              "pryr", "dplyr", "ggplot2", "cowplot",
              "pheatmap", "igraph", "RColorBrewer",
              "BiocParallel", "DESeq2", "limma", 
              "org.Mm.eg.db", "scales",
              "gatom", "rUtils", "GAM", 
              "GAM.db", "GAM.networks")
lapply(Packages, library, character.only = TRUE)

# devtools::install_github("ctlab/gatom")
# devtools::install_github("ctlab/GAM")
# devtools::install_github("assaron/rUtils")
# source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")
# load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))

data("met.id.map")
data("kegg.mouse.network")
GAM:::lazyData("met.id.map")
GAM:::lazyData("kegg.db")
enz2gene <- kegg.db$enz2gene
rxn2enz <- kegg.db$rxn2enz
load("Data/reflink.rda")

tissueColorsList <- list(
  Aorta = "#771155",
  Blood = "#AA4488",
  `Bone marrow` = "#CC99BB",
  Brain = "#114477",
  `Embryonic body` = "hotpink",
  Colon = "#77AADD", # "#4477AA"
  Dermis = "#117744",
  `Dorsal root ganglia` = "#44AA77",
  `Epididymal fat` = "#88CCAA",
  Epithelium = "#777711",
  Heart = "#AAAA44",
  `Inguinal fat` = "#DDDD77",
  Kidney = "#774411",
  Liver = "lightseagreen",
  Lung = "darkviolet", # "#AA7744"
  `Mediastinal LN` = "#771122",
  `Mesenteric fat` = "#AA4455",
  `Mesenteric LN` = "#DD7788",
  `Mesenteric sheet` = "gray0",
  `Peritoneal cavity` = "gray40",
  `Peyer's patch` = "lawngreen",
  `Sciatic nerve` = "yellow",
  `Spinal cord` = "cyan",
  Spleen = "darkgoldenrod1",
  Thymus = "darkgoldenrod4",
  `Yolk sac` = "red")

metaColorsList <- list(
  `tissue MF` = "yellow4",
  EMP = "hotpink",
  `MF (E6-E8)` = "steelblue1",
  `tissue DC` = "red",
  Mo = "lightpink4",
  pDC = "green3",
  MG = "tan2",
  migDC = "purple",
  `alv MF & SPM` = "forestgreen",
  `adipose MF` = "darkslategray")

cellColorsList <- list(
  `Dendritic cell` = "sandybrown",
  `Erythro-myeloid progenitor` = "red",
  `Macrophage (E6-E8)` = "deeppink2",
  Microglia = "turquoise",
  Macrophage = "purple",
  Monocyte = "forestgreen")

cellColorsList2 <- list(
  cDC1 = "darkgoldenrod4",
  cDC2 = "yellow3",
  DC = "sandybrown",
  EMP = "red",
  MF = "purple",
  `MF (E14.5-E18.5)` = "gray1",
  `MF (E6-E8)` = "deeppink2",
  MG = "turquoise",
  Mo = "forestgreen",
  pDC = "brown")

# shapes <- c("circle plus",
#             "circle",
#             "circle cross",
#             "cross",
#             "asterisk",
#             "triangle",
#             "triangle down filled",
#             "star",
#             "square",
#             "square triangle",
#             "square cross",
#             "square plus",
#             "plus",
#             "diamond filled",
#             "diamond plus",
#             "diamond filled")

pcaPlot <- function (es, c1, c2, scale = F, size = 2)
{
  stopifnot(require(ggplot2))
  pca <- prcomp(t(exprs(es)), scale. = scale)
  explained <- (pca$sdev)^2/sum(pca$sdev^2)
  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)
  a <- cbind(as.data.frame(pca$x), the_sample = colnames(es),
             pData(es))
  pp <- ggplot(data = a)
  pp + aes_string(x = xs[c1], y = xs[c2]) + geom_point(size = size) +
    xlab(xlabs[c1]) + ylab(xlabs[c2])
}

gradientColor <- c("blue", "grey75", "red")

"%o%" <- compose

scriptName <- function() sub(".*=", "", head(grep("^--file=", commandArgs(), value=T), 1))

if (!interactive()) {
  options(error = quote({
    saveTo <- paste0(scriptName(), ".dump")
    dump.frames(saveTo,
                to.file=TRUE,
                include.GlobalEnv = TRUE)
    quit(save="no", status=1)
  }))

  mypdf <- function(filename = paste0(scriptName(), ".pdf"), ...) {
    pdf(filename, width=8, height=6, ...)
  }
  options(device="mypdf")
}

processModule <- function(module, dir, s=34, name=NULL) {
  if (is.null(name)) {
    name <- deparse(substitute(module))
  }
  #     module <- expandReactionNodeAttributesToEdges(module)
  #     E(module)$label <- ""
  if (is.null(V(module)$logPval)) {
    V(module)$logPval <- -5
  }

  if (is.null(V(module)$log2FC)) {
    V(module)$log2FC <- NA
  }

  dir.create(dir, showWarnings=F)

  file <- file.path(dir, paste0(name, ".dot"))
  saveModuleToDot(module, file=file, name=name)
  saveModuleToXgmml(module, name=name, file=file.path(dir, paste0(name, ".xgmml")))
  system2("neato", c("-Tpdf",
                     "-o", file.path(dir, paste0(name, ".pdf")),
                     file), stderr = F)
  system2("neato", c("-Tpng",
                     "-o", file.path(dir, paste0(name, ".png")),
                     file), stderr = F)
  gatom::saveModuleToPdf(module = module, seed = s, file = paste0(dir, name, "_nice", s, ".pdf"))
  }

makeAtomGraph <- function(network,
                          org.gatom.anno,
                          gene.de,
                          gene.de.meta=getGeneDEMeta(gene.de, org.gatom.anno),
                          gene.keep.top=12000,
                          met.db,
                          met.de,
                          met.de.meta=getMetDEMeta(met.de, met.db)) {
  if (!is.null(gene.de)) {
    gene.de <- prepareDE(gene.de, gene.de.meta)
    gene.de <- gene.de[signalRank <= gene.keep.top]
  }

  met.de <- prepareDE(met.de, met.de.meta)

  edge.table <- .makeEdgeTable(network=network,
                               org.gatom.anno=org.gatom.anno,
                               gene.de=gene.de,
                               gene.de.meta=gene.de.meta)
  all.atoms <- union(edge.table$atom.x, edge.table$atom.y)
  vertex.table <- .makeVertexTable(network=network,
                                   atoms=all.atoms,
                                   met.db=met.db,
                                   met.de=met.de,
                                   met.de.meta=met.de.meta)
  g <- graph.data.frame(edge.table, directed=FALSE, vertices = vertex.table)
  gc <- components(g)
  g <- induced.subgraph(g, gc$membership == which.max(gc$csize))

  if (!is.null(met.to.filter)) {
    nodes.to.del <- V(g)[metabolite %in% met.to.filter]
    if (length(nodes.to.del) == 0) {
      warning("Found no metabolites to mask")
    } else {
      g <- delete_vertices(g, v = V(g)[metabolite %in% met.to.filter])
    }
  }
  g
}

makeENetworkWithScore <- function(score, base, network) {
  fake.gene.de <- data.frame(
    ID=names(score),
    pval=2 ^ -(score),
    score=score,
    log2FC=score)
  fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]),
                        by.x="ID", by.y="Entrez", all.x=T)

  suppressWarnings(es.re <- makeExperimentSet(network, gene.de=fake.gene.de,
                                              reactions.as.edges=T, plot=F))
  subnet.scored <- es.re$subnet

  #V(subnet.scored)$score <- -0.1
  V(subnet.scored)$score <- -0.01
  rxn.vertices <- V(subnet.scored)[nodeType == "rxn"]
  #     E(subnet.scored)$score <-
  #         -(1 - pmax(1 - E(subnet.scored)$dist.closest, 0)) +
  #         base
  #

  es.re.scored <- es.re
  es.re.scored$subnet.scored <- subnet.scored
  es.re.scored
}

sgmwcs.batchSolver2 <- function(sgmwcs, nthreads = 1, timeLimit = -1,
                                nodes.group.by=NULL, edges.group.by=NULL,  group.only.positive=F,
                                c.size=NULL,
                                quite=T) {
  function(networks) {
    N <- length(networks)
    graphBatch.dir <- tempfile("graphBatch")
    dir.create(graphBatch.dir)
    graph.dirs <- file.path(graphBatch.dir, seq_len(N))

    instances <- lapply(seq_len(N), function(i) {
      network <- networks[[i]]
      score.edges <- "score" %in% list.edge.attributes(network)
      score.nodes <- "score" %in% list.vertex.attributes(network)
      if (!score.nodes) {
        V(network)$score <- 0
      }
      if (!score.edges) {
        E(network)$score <- 0
      }

      graph.dir <- graph.dirs[i]
      dir.create(graph.dir)

      gatom:::writeSgmwcsInstance(graph.dir = graph.dir,
                                  network=network,
                                  nodes.group.by = nodes.group.by,
                                  edges.group.by = edges.group.by,
                                  group.only.positive = group.only.positive
      )
    })

    system2(sgmwcs, c("-a", sprintf("1-%s", N),
                      "-d", graphBatch.dir,
                      "--threads", nthreads,
                      "--timelimit", timeLimit,
                      if (!is.null(c.size)) c("-c", c.size) else NULL,
                      "--break"
    ))

    res <- lapply(seq_len(N), function(i) {
      network <- networks[[i]]
      instance <- instances[[i]]
      graph.dir <- graph.dirs[i]
      solution.file <- paste0(instance$nodes.file, ".out")

      if (!file.exists(solution.file)) {
        warningf("Solution file '%s' not found", solution.file)
        NULL
      } else {
        gatom:::readGraph(
          node.file = solution.file,
          edge.file = paste0(instance$edges.file, ".out"),
          network = network)
      }
    })
  }
}

makeNetwork <- function(gene.de, met.de=NULL, largest.component=TRUE) {    
    stopifnot(require(igraph))
    GAM:::lazyData("kegg.db")
    
    if (! "logPval" %in% names(gene.de)) {
        gene.de$logPval <- log2(gene.de$pval)    
    }
    
    if (!is.null(met.de) && ! "logPval" %in% names(met.de)) {
        met.de$logPval <- log2(met.de$pval)    
    }
    
    enz.de <- convertPval(gene.de, kegg.db$enz2gene$gene, kegg.db$enz2gene$enz)
    rxn.de <- convertPval(enz.de, kegg.db$rxn2enz$enz, kegg.db$rxn2enz$rxn)
    rpair.de <- convertPval(rxn.de, rpairs$rxn, rpairs$rpair)
    
    
    edges <- merge(rpaligns, rpair.de, by.x="rpair", by.y="ID")
    #edges
    #edges[atom.x > atom.y,] list(atom.x, atom.y)]  <- edges[atom.x > atom.y, list(atom.y, atom.x)]
    edges <- edges[!duplicated(edges[, c("atom.x", "atom.y")]), ]
    edges <- edges[edges$atom.x != edges$atom.y, ]
    edges$label <- edges$symbol
    
    nodes <- data.frame(ID=with(rpaligns, unique(c(atom.x, atom.y))), stringsAsFactors=F)
    nodes$KEGG <- gsub("_.*$", "", nodes$ID)
    nodes <- merge(nodes, kegg2name, by.x="KEGG", by.y="met", all.x=T)
    if (!is.null(met.de)) {
        nodes <- merge(nodes, met.de, by.x="KEGG", by.y="ID", all.x=T, suffixes=c("", ".orig"))
    }
    #nodes <- rename(nodes, c("name"="label"))
    
    
    atom.graph <- GAM:::graph.from.tables(node.table=nodes, node.col="ID", name.as.label=T,
                                          edge.table=edges, 
                                          edge.cols=c("atom.x", "atom.y"), directed=F)
    
    
    res <- delete.vertices(atom.graph, V(atom.graph)[KEGG %in% kegg.db$mets2mask]$name)
    if (largest.component) {
        gc <- clusters(res)
        res <- induced.subgraph(res, gc$membership == which.max(gc$csize))

    }  
        
    res
}

collapseAtoms2 <- function(m) {
    t <- data.frame(v=seq_along(V(m)), met=gsub("_.*$", "", V(m)$name))
    t <- data.frame(v=V(m)$name, met=gsub("_.*$", "", V(m)$name), stringsAsFactors=F)
    toCollapse <- merge(t, t, by="met")
    toCollapse <- toCollapse[(toCollapse$v.x < toCollapse$v.y), ]
    
    res <- igraph::add.edges(m, matrix(c(toCollapse$v.x, toCollapse$v.y), nrow=2, byrow=T))    
    res
}

collapseAtoms <- function(m) {
    nodes <- data.frame(GAM:::get.vertex.attributes(m), stringsAsFactors=F)
    nodes$name <- gsub("_.*$", "", nodes$name)
    nodes <- unique(nodes)
    edges <- cbind(get.edgelist(m), data.frame(GAM:::get.edge.attributes(m)))
    edges[, c(1,2)] <- lapply(edges[, 1:2], gsub, pattern="_.*$", replacement="")
    edges <- unique(edges)
    res <- GAM:::graph.from.tables(node.table=nodes, edge.table=edges, name.as.label=F)
    res
}

replaceNA <- function(x, y) {
    ifelse(is.na(x), y, x)
}

addEdges <- function(m, net, es) {    
    m.temp <- induced.subgraph(net, V(m)$name)
    
    m.edges <- get.edge.attributes(m)
    m.temp.edges <- get.edge.attributes(m.temp)
    
    edges.old.id <- which(duplicated(rbind(m.edges, m.temp.edges))) - nrow(m.edges)
        
    stopifnot(length(edges.old.id) == length(E(m)))
            
    m.new <- subgraph.edges(m.temp, unique(c(edges.old.id, E(m.temp)[eval(es)])))
    m.new
}

processModule <- function(module, dir, name=NULL, do.png=F, keep.dot=F) {
    if (is.null(name)) {
        name <- deparse(substitute(module))
    }
#     module <- expandReactionNodeAttributesToEdges(module)
#     E(module)$label <- ""
    if (is.null(V(module)$logPval)) {
        V(module)$logPval <- -5    
    }
    
    if (is.null(V(module)$log2FC)) {
        V(module)$log2FC <- NA   
    }
    
    dir.create(dir, showWarnings=F)
    
    file <- file.path(dir, paste0(name, ".dot"))
    saveModuleToDot(module, file=file, name=name)
    saveModuleToXgmml(module, name=name, file=file.path(dir, paste0(name, ".xgmml")))
    system2("neato", c("-Tpdf", 
                       "-o", file.path(dir, paste0(name, ".pdf")),
                       file), stderr = F)
    if (do.png) {
        system2("neato", c("-Tpng", 
                           "-o", file.path(dir, paste0(name, ".png")),
                           file), stderr = F)    
    }
    if (!keep.dot) {
        file.remove(file)
    }
    
}

mwcsize <- function(g) {         
    was.connected <- is.connected(g)
    et <- as.data.table(get.edge.attributes(g))
    et <- et[, list(from, to, score)]
    et[, from.weight := pmax(-V(g)[et$from]$score, 0)+1e-3]
    et[, to.weight := pmax(-V(g)[et$to]$score, 0)+1e-3]
    et <- et[score > 0, ]
    et[, from.d := score * from.weight / (from.weight + to.weight)]
    et[, to.d := score * to.weight / (from.weight + to.weight)]
    
    ds <- rbind(
        et[, list(v=from, d=from.d)],
        et[, list(v=to, d=to.d)])
    
    ds <- aggregate(d ~ v, data=ds, sum)
    
    V(g)[ds$v]$score <- V(g)[ds$v]$score + ds$d
    E(g)$origEdge <- E(g)
    E(g)[score > 0]$score <- 0
    
    
    neg.edges <- E(g)[score < 0]
    neg.edges.table <- get.edgelist(g)[neg.edges, ]
    neg.edges.names <- sprintf("%s_%s_%s", neg.edges.table[,1], neg.edges.table[,2], neg.edges)
    V(g)$wasEdge <- FALSE
    g <- add.vertices(g, length(neg.edges), name=neg.edges.names, score=E(g)[neg.edges]$score, wasEdge=TRUE, origEdge=neg.edges)
    new.edges <- rbind(
        cbind(neg.edges.table[,1], neg.edges.names),
        cbind(neg.edges.table[,2], neg.edges.names)
    )
    g <- add.edges(g, t(as.matrix(new.edges)), score=0)
    g <- delete.edges(g, neg.edges)    
    E(g)$score <- 0
    stopifnot(is.connected(g) == was.connected)
    g
}

demwcsize <- function(m, g) {
    E(g)$origEdge <- E(g)
    orig.vertices <- V(g)[name %in% V(m)$name]
    edges.to.add <- unique(c(
        E(induced.subgraph(g, orig.vertices))[score >= 0]$origEdge,
        na.omit(c(V(m)$origEdge, E(m)$origEdge))))        
    
    m.x <- subgraph.edges(g, edges.to.add)    
    m.x <- remove.edge.attribute(m.x, "origEdge")
    m.x
}

get.score <- function(m) sum(V(m)$score) + sum(E(m)$score)

heinz2mc.solver <- function (heinz2, niter = 1000) 
{
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (score.edges) {
            network <- GAM:::MWCSize(network.orig)
        }
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        solution.file <- file.path(graph.dir, "sol.txt")
        system2(paste0(heinz2), c("-n", nodes.file, "-e", edges.file, 
                                  "-o", solution.file, "-m", niter, "-v", 
                                  1))
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- BioNet::readHeinzGraph(node.file = solution.file, network = network, 
                              format = "igraph")
        if (score.edges) {
            res <- GAM:::deMWCSize(res, network.orig)
        }
        return(res)
    }
}

heinz2.solver <- function (heinz2, nthreads = 1, timeLimit = -1) 
{
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (score.edges) {
            network <- GAM:::MWCSize(network.orig)
        }
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        solution.file <- file.path(graph.dir, "sol.txt")
        system2(paste0(heinz2), c("-n", nodes.file, "-e", edges.file, 
                                  "-o", solution.file, "-m", nthreads, "-v", 
                                  0, "-t", timeLimit))
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- BioNet::readHeinzGraph(node.file = solution.file, network = network, 
                                      format = "igraph")
        if (score.edges) {
            res <- GAM:::deMWCSize(res, network.orig)
        }
        return(res)
    }
}

writeGmwcsInstance <- function(graph.dir, network,
                   nodes.group.by=NULL, 
                   edges.group.by=NULL,
                   group.only.positive=F) {
    
    dir.create(graph.dir, showWarnings = FALSE)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    synonyms.file <- file.path(graph.dir, "synonyms.txt")
    
    synonyms <- c()
    
    nt <- get.vertex.attributes(network)
    if (!is.null(nodes.group.by)) {
        f <- as.formula(sprintf("name ~ %s", nodes.group.by))
        if (all.vars(f) %in% colnames(nt)) {
            synonyms <- c(synonyms, aggregate(f, data=nt, paste0, collapse=" ")$name)    
        } else {
            warningf("Can't collapse nodes, not all fields present: %s", 
                     paste0(setdiff(all.vars(f), colnames(nt)), collapse=", "))
            synonyms <- c(synonyms, nt$name)    
        }
        
    } else {
        synonyms <- c(synonyms, nt$name)
    }
    nt <- rename(nt[, c("name", "score")], c("name"="#name"))
    
    et <- get.edge.attributes(network, include.ends = T)
    if (!is.null(edges.group.by)) {
        etx <- if (group.only.positive) { 
            synonyms <- c(synonyms, with(et[et$score <= 0,], sprintf("%s -- %s", from, to)))
            et[et$score > 0,] 
        } else { 
            et 
        }
        if (nrow(etx) > 0) {
            synonyms <- c(synonyms, 
                          aggregate(name ~ edges.group.by,
                                    data=list(
                                        name=sprintf("%s -- %s", etx$from, etx$to),
                                        edges.group.by=etx[[edges.group.by]]),
                                    paste0, collapse=" ")$name
            )    
        }
    }        
    et <- rename(et[, c("from", "to", "score")], c("from"="#from"))
    
    write.tsv(nt, file=nodes.file)
    write.tsv(et, file=edges.file)
    writeLines(sprintf("%s", synonyms), con=synonyms.file)
    
    if (length(synonyms) == 0) {
        synonyms.file <- NULL
    }
    
    list(nodes.file=nodes.file, 
         edges.file=edges.file, 
         synonyms.file=synonyms.file)
}

gmwcs.solver <- function (gmwcs, nthreads = 1, timeLimit = -1, nodes.group.by=NULL, edges.group.by=NULL, c.size=NULL, group.only.positive=F, quite=T,
                          minimize.size=F) {
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (!score.edges) {
            E(network)$score <- 0
        }
        
        graph.dir <- tempfile("graph")
        
        instance <- writeGmwcsInstance(graph.dir = graph.dir,
                                       network=network,
                                       nodes.group.by = nodes.group.by,
                                       edges.group.by = edges.group.by,
                                       group.only.positive = group.only.positive
        )
        
        system2(gmwcs, c("--nodes", instance$nodes.file,
                                "--edges", instance$edges.file,
                                if (!is.null(instance$synonyms.file)) c("--synonyms", instance$synonyms.file, "-B", 1) else NULL,
                                "--threads", nthreads, 
                                "--timelimit", timeLimit,
                                if (!is.null(c.size)) c("-c", c.size) else NULL,
                                if (minimize.size) c("-p", 1e-3) else NULL,
                                "--break"
        ))
        solution.file <- paste0(instance$nodes.file, ".out")
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- GAM:::readGraph(node.file = solution.file,
                               edge.file = paste0(instance$edges.file, ".out"),
                               network = network)
        #attr(res, "optimal") <- any(grepl("SOLVED TO OPTIMALITY", out))
        return(res)
    }
}

writeGraph2 <- function(g, prefix) {
    nt <- get.vertex.attributes(g)
    et <- get.edge.attributes(g)
    nt <- nt[, c("name", "score", "KEGG", "formula")]
    et <- et[, c("from", "to", "score", "rpair", "symbol")]    
    write.tsv(nt, file=paste0(prefix, ".nodes.tsv"))
    write.tsv(et, file=paste0(prefix, ".edges.tsv"))
}

doGatom <- function(dir, tag, gene.de, met.de, k=100, 
                    e.threshold.min=-log2(0.01),
                    v.threshold.min=-log2(0.01),
                    del=NULL) {
    dir.create(dir, showWarnings = F)
    
    if (length(del) > 0) {
        tag <- paste0(tag, paste0("_d", del, collapse=""))
    }
    
    a <- makeNetwork(gene.de=gene.de, met.de=met.de)
    g1 <- a
    
    if (!is.null(met.de)) {
        v.threshold <- sort(unique(V(g1)$logPval))[k]
        # v.threshold <- -log2(0.01)
        v.threshold <- max(v.threshold, v.threshold.min)
        messagef("Metabolite threshold: %f", v.threshold)    
        V(g1)$score <- with(get.vertex.attributes(g1), -replaceNA(logPval, 0) - v.threshold)
    } else {
        V(g1)$score <- 0
    }

    if (!"score" %in% list.edge.attributes(g1)) {
        # e.threshold <- -log2(0.01)
        e.threshold <- -sort(unique(E(g1)$logPval))[k]
        e.threshold <- max(e.threshold, e.threshold.min)
        messagef("Gene threshold: %f", e.threshold)
        E(g1)$score <- with(get.edge.attributes(g1), -logPval - e.threshold)
        E(g1)$score <- pmin(E(g1)$score, 1000)
    } else {
        messagef("Genes with positive score: %s", sum(unique(E(g1)$score) > 0))
    }
    
    g1 <- delete.edges(g1, E(g1)[origin %in% del])
    
    de.mod <- solve(g1)
    
    de.mod1 <- de.mod
    
    if (is.null(met.de)) {
        V(de.mod1)$score <- -0.01
    }
  
    de.mod2 <- solve.syn2(de.mod1)
    
    processModule(collapseAtoms2(de.mod2), dir=dir, 
                  name=paste0(tag, 
                              ".k", k,
                              if (!is.null(gene.de)) { ".g" } else NULL,  
                              if (!is.null(met.de)) { ".m" } else NULL,  
                              ".col2"))
    
    processModule(collapseAtoms(de.mod2), dir=dir, 
                  name=paste0(tag, 
                              ".k", k,
                              if (!is.null(gene.de)) { ".g" } else NULL,  
                              if (!is.null(met.de)) { ".m" } else NULL,  
                              ".col"))
    
    et <- unique(get.edge.attributes(g1, include.ends = T)[, c("symbol", "pval", "score")])
    et <- et[order(et$pval), ]
    et10 <- head(et, n=10)
    print(et10[!et10$symbol %in% E(de.mod2)$symbol, ])
    
    de.mod2
}

doGatom2 <- function(solve, dir, tag, gene.de, met.de, k=100, k.met=k, 
                    e.threshold.min=0.1,
                    v.threshold.min=0.1,
                    del=NULL,
                    add.top=3000,
                    remove.zeros=T,
                    do.col2=F,
                    do.png=F) {
    if (!is.null(dir)) {
        dir.create(dir, showWarnings = F)
    }
    
    if (length(del) > 0) {
        tag <- paste0(tag, paste0("_d", del, collapse=""))
    }
    
    a <- makeNetwork(gene.de=gene.de, met.de=met.de)
    g1 <- a
    
    if (!is.null(met.de)) {
        #met.bum <- fitBumModel(unique(met.de$pval), plot = F)
        met.bum <- fitBumModel(na.omit(unique(V(g1)$pval)), plot = F)
        v.threshold <- sort(unique(V(g1)$pval))[k.met]
        v.threshold <- min(v.threshold, fdrThreshold(v.threshold.min, met.bum))
        messagef("Metabolite threshold: %f", v.threshold)    
        V(g1)$score <- with(get.vertex.attributes(g1), 
                            (met.bum$a - 1) *
                            (log(replaceNA(pval, 1)) - log(v.threshold)))
    } else {
        V(g1)$score <- 0
    }
    
    if (!"score" %in% list.edge.attributes(g1)) {
        #gene.bum <- fitBumModel(unique(gene.de$pval), plot = F)
        gene.bum <- fitBumModel(na.omit(unique(E(g1)$pval)), plot = F)
        
        if (gene.bum$a > 0.6) {
            messagef("Gene BUM alpha is too big: %f", gene.bum$a)    
            return(NULL)
        }
        
        # e.threshold <- -log2(0.01)
        e.threshold <- sort(unique(E(g1)$pval))[k]
        e.threshold <- min(e.threshold, fdrThreshold(e.threshold.min, gene.bum))
        
        messagef("Gene threshold: %f", e.threshold)
        E(g1)$score <- with(get.edge.attributes(g1), 
                            (gene.bum$a - 1) *
                                (log(replaceNA(pval, 1)) - log(e.threshold)))
        
        E(g1)$score <- pmin(E(g1)$score, 1000)
    } else {
        messagef("Genes with positive score: %s", sum(unique(E(g1)$score) > 0))
    }
    
    if (!is.null(del)) {
        g1 <- delete.edges(g1, E(g1)[origin %in% del])    
    }
    
    de.mod <- solve(g1)
    
    de.mod1 <- de.mod
    
    if (remove.zeros) {
        de.mod2 <- solve(de.mod1, minimize.size=TRUE)
    } else {
        de.mod2 <- de.mod1
    }
    
    if (!is.null(gene.de) && add.top > 0) {
        superExpressed <- gene.de[head(order(gene.de$baseMean, decreasing = T), n=add.top), ]$ID
        
        toAdd <- data.table(get.edge.attributes(g1, index = E(g1)[origin %in% superExpressed], include.ends = T))
        
        toAdd <- toAdd[from %in% V(de.mod2)$name & to %in% V(de.mod2)$name]
        
        toAdd <- toAdd[!paste0(from, ".", to) %in% with(get.edge.attributes(de.mod2, include.ends = T), paste0(from, ".", to))]
        
        de.mod2 <- add.edges(de.mod2, rbind(toAdd$from, toAdd$to), attr=tail(as.list(toAdd), -2))
    }
    
    if (!is.null(dir)) {
        if (do.col2) {
            processModule(collapseAtoms2(de.mod2), dir=dir, 
                          name=paste0(tag, 
                                      ".k", k,
                                      if (!is.null(gene.de)) { ".g" } else NULL,  
                                      if (!is.null(met.de)) { ".m" } else NULL,  
                                      ".col2"))
        }
        
        processModule(collapseAtoms(de.mod2), dir=dir, 
                      name=paste0(tag, 
                                  ".k", k,
                                  if (!is.null(gene.de)) { ".g" } else NULL,  
                                  if (!is.null(met.de)) { ".m" } else NULL,  
                                  ".col"))
    }
    
    et <- unique(get.edge.attributes(g1, include.ends = T)[, c("symbol", "pval", "score")])
    et <- et[order(et$pval), ]
    et10 <- head(et, n=10)
    print(et10[!et10$symbol %in% E(de.mod2)$symbol, ])
    
    de.mod2
}

gseaFisher <- function(pathways, universe, genes, minSize=10, maxSize=Inf) {
    pathways <- lapply(pathways, intersect, universe)
    pathways <- pathways[sapply(pathways, length) >= minSize & sapply(pathways, length) <= maxSize]
    genes <- intersect(genes, universe)    
    
    overlap <- lapply(lapply(pathways, intersect, genes), sort)
    overlaps <- data.table(
        q=sapply(overlap, length),
        m=sapply(pathways, length),
        n=length(universe)-sapply(pathways, length),
        k=length(genes))
    
    # q-1 because we want probability of having >=q white balls
    pathways.pvals <- with(overlaps,
                           mapply(phyper, q-1, m, n, k, lower.tail = FALSE))
    res <- data.table(pathway=names(pathways), 
                      pval=pathways.pvals,
                      k=overlaps$q,
                      K=overlaps$m)
    res[, padj := p.adjust(pval, method="BH")]
    res[, overlap := .(overlap)]    
    
    
    res <- res[order(pval),]
    res[]
    res
}

heinz2mc.solver <- function (heinz2, niter = 1000)
{
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (score.edges) {
            network <- GAM:::MWCSize(network.orig)
        }
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        solution.file <- file.path(graph.dir, "sol.txt")
        system2(paste0(heinz2), c("-n", nodes.file, "-e", edges.file,
                                  "-o", solution.file, "-m", niter, "-v",
                                  0))
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- BioNet::readHeinzGraph(node.file = solution.file, network = network,
                                      format = "igraph")
        if (score.edges) {
            res <- GAM:::deMWCSize(res, network.orig)
        }
        return(res)
    }
}

heinz2.solver <- function (heinz2, nthreads = 1, timeLimit = -1)
{
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (score.edges) {
            network <- GAM:::MWCSize(network.orig)
        }
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        solution.file <- file.path(graph.dir, "sol.txt")
        system2(paste0(heinz2), c("-n", nodes.file, "-e", edges.file,
                                  "-o", solution.file, "-m", nthreads,
                                  "-v", 0, "-t", timeLimit,
                                  "-p"
                                  ))
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- BioNet::readHeinzGraph(node.file = solution.file, network = network,
                              format = "igraph")
        if (score.edges) {
            res <- GAM:::deMWCSize(res, network.orig)
        }
        return(res)
    }
}

get.score <- function(m) sum(V(m)$score) + sum(E(m)$score)

makeNetworkWithDist <- function(gene.dist.to.cluster, base, network) {
    fake.gene.de <- data.frame(
        ID=names(gene.dist.to.cluster),
        pval=2 ^ ((gene.dist.to.cluster-2)*20),
        dist.closest=gene.dist.to.cluster,
        log2FC=(1-gene.dist.to.cluster)*2)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    es.rn <- makeExperimentSet(network, gene.de=fake.gene.de, reactions.as.edges=F, plot=F)
    subnet.scored <- es.rn$subnet

    #V(subnet.scored)$score <- -0.1
    V(subnet.scored)$score <- -0.01
    rxn.vertices <- V(subnet.scored)[nodeType == "rxn"]
    V(subnet.scored)[rxn.vertices]$score <-
        -(1 - pmax(1 - V(subnet.scored)[rxn.vertices]$dist.closest, 0)) +
        base

    repeats <- table(V(subnet.scored)[score > 0]$origin)
    repeats <- names(repeats[repeats >= 4])
    V(subnet.scored)[origin %in% repeats]$score <- -0.01


    es.rn.scored <- es.rn
    es.rn.scored$subnet.scored <- subnet.scored
    es.rn.scored
}

makeENetworkWithDist <- function(gene.dist.to.cluster, base, network) {
    fake.gene.de <- data.frame(
        ID=names(gene.dist.to.cluster),
        pval=2 ^ ((gene.dist.to.cluster-2)*20),
        dist.closest=gene.dist.to.cluster,
        log2FC=(1-gene.dist.to.cluster)*2)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    es.re <- makeExperimentSet(network, gene.de=fake.gene.de, reactions.as.edges=T, plot=F)
    subnet.scored <- es.re$subnet

    #V(subnet.scored)$score <- -0.1
    V(subnet.scored)$score <- -0.01
    rxn.vertices <- V(subnet.scored)[nodeType == "rxn"]
    E(subnet.scored)$score <-
        -(1 - pmax(1 - E(subnet.scored)$dist.closest, 0)) +
        base


    es.re.scored <- es.re
    es.re.scored$subnet.scored <- subnet.scored
    es.re.scored
}

getCenter <- function(gene.exprs, cluster.genes=seq_len(nrow(gene.exprs)), cluster.genes.neg=c(), method=c("pearson", "spearman")) {
    method <- match.arg(method)

    cluster.exprs <- rbind(
        gene.exprs[cluster.genes,, drop=F],
        -gene.exprs[cluster.genes.neg,, drop=F])

    if (method == "spearman") {
        cluster.exprs <- t(apply(cluster.exprs, 1, rank))
    }

    cluster.exprs.znorm <- zScore(cluster.exprs)

    center <- apply(cluster.exprs.znorm, 2, mean, na.rm=T)
    center <- center / sd(center, na.rm = T)
    center[!is.finite(center)] <- NA
    center
}

gamCluster <- function(gene.exprs, tag, work.dir, network=kegg.mouse.network, show_colnames=T) {
    stopifnot(require(cluster))
    stopifnot(require(GAM))
    stopifnot(require(igraph))
    stopifnot(require(pheatmap))
    stopifnot(require(data.table))
    dir.create(work.dir, showWarnings=F, recursive=T)
    gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$gene.id.map$Entrez, ]
    gene.exprs <- gene.exprs[apply(gene.exprs, 1, sd) > 0,]

    gene.cor <- cor(t(gene.exprs))
    gene.cor.dist <- as.dist(1 - gene.cor)
    gene.cor.dist.m <- as.matrix(gene.cor.dist)

    gK <- 2
    gK <- 10
    reorder <- sample(gK)
    set.seed(42)
    gene.pam <- pam(gene.cor.dist, k=gK)
    if (FALSE) {
        pheatmap(
            normalize.rows(gene.exprs[order(reorder[gene.pam$clustering]), ]),
            cluster_rows=F, cluster_cols=T,
            show_rownames=F, show_colnames=show_colnames)
        pheatmap(
            normalize.rows(gene.exprs[gene.pam$medoids, ]),
            cluster_rows=F, cluster_cols=T,
            show_rownames=F, show_colnames=show_colnames)

    }

    gene.dist.to.clusters <- as.matrix(gene.cor.dist)[, gene.pam$medoids]
    gene.dist.to.cluster.closest <- apply(gene.dist.to.clusters, 1, min)
    gene.dist.to.cluster.which <- apply(gene.dist.to.clusters, 1, which.min)


    modules.gene <- list()
    modules.gene.centers <- list()

    gene.bases <- c(rep(0.025, gK))
    i <- 10

    #dir.create(sprintf("%s/refined", work.dir), showWarnings=FALSE, recursive=TRUE)

    for (i in seq_len(gK)) {
        cluster.genes <- names(gene.pam$clustering)[gene.pam$clustering == i]
        #weights <- rep(1, length(cluster.genes))

        medoid <- gene.exprs[gene.pam$medoids[i],]
        medoid <- (medoid - mean(medoid)) / sd(medoid)

        it = 0
        base <- gene.bases[i]


        centers <- rbind(medoid)
        rownames(centers)[1] <- paste0("medoid.", i)

        redo <- TRUE
        while (TRUE) {
            center <- getCenter(gene.exprs, cluster.genes)

            cor.diff <- 1 - max(apply(centers, 1, cor, y=center))
            difference <- min(apply(sweep(centers, 2, center), 1, function(x) { sum(abs(x)) }))
            centers <- rbind(center, centers)
            rownames(centers)[1] <- paste0("center.", i, ".r", it)
            if (FALSE) {
                pheatmap(
                    centers,
                    cluster_rows=F, cluster_cols=F,
                    show_colnames=T)
            }

            if (!redo && cor.diff < 0.001) {
                break
            }
            redo <- FALSE

            print(paste0("Difference: ", cor.diff))
            print(paste0("Core: ", paste0(core, collapse=" ")))
            message(sprintf("iteration #%s, base=%s", it, base))

            gene.dist.to.cluster <- 1 - cor(t(gene.exprs), as.matrix(center))[,1]
            es.rn.scored <- makeNetworkWithDist(gene.dist.to.cluster, base, network)


            #module <- findModule(es.rn.scored, solver.mc)
            module <- findModule(es.rn.scored, solver.heinz2)

            module

            if (is.null(module)) {
                module <- findModule(es.rn.scored, solver.heinz)
            }

            if (FALSE) {
                n <- es.rn.scored$subnet.scored
                V(n)$inModule <- F
                V(n)[name %in% V(module)$name]$inModule <- T
                vInModule <- V(n)[inModule]
                E(n)$inModule <- F
                E(n)[vInModule %--% vInModule]$inModule <- T
                tag <- paste0("sol.", i, ".r", it, ".b", base)
                saveModuleToXgmml(n, file=file.path(work.dir, paste0(tag, ".xml")), name=tag)
            }

            modules.gene[[i]] <- module
            modules.gene.centers[[i]] <- centers

            #V(module)$betweenness <- betweenness(module)
            #bet.median <- median(V(module)$betweenness)
            #core <- unique(na.omit(V(module)[betweenness < bet.median]$origin))
            core <- unique(na.omit(V(module)$origin))
            cluster.genes.new <- intersect(core, rownames(gene.exprs))
            if (length(cluster.genes.new) >= 10) {
                cluster.genes <- cluster.genes.new
            }
            
            #         cluster.genes<- intersect(E(module)$RefSeq, rownames(gene.exprs))
#             saveModuleToXgmml(
#                 module,
#                 file=sprintf("%s/refined/c%s.gene.%s.r%s.xgmml", work.dir, gK, i, it),
#                 name=sprintf("%s.gene.%s.refined.r%s.b%s", tag, i, it, base))
            if (length(V(module)) < 50 || get.score(module) <= 0 || length(core) < 10) {
                base <- base + 0.025
                redo <- TRUE
            } else if (length(V(module)) > 150) {
                base <- base - 0.05
                redo <- TRUE
            }
            it = it + 1
        }

#         pheatmap(
#             centers,
#             cluster_rows=F, cluster_cols=T,
#             show_colnames=F)
#
        saveModuleToXgmml(
            module,
            file=sprintf("%s/c%s.gene.%s.xgmml", work.dir, gK, i),
            name=sprintf("%s.gene.%s.refined.b%s", tag, i, base))
    }

#     pheatmap(
#         do.call(rbind, modules.gene.centers),
#         cluster_rows=F, cluster_cols=T,
#         file=sprintf("%s/centroids.gene.all.png", work.dir))
#
    res <- list(modules=modules.gene, centers=modules.gene.centers)
    writeClusterGenes(exprs=gene.exprs, cluster=res, tag=tag, work.dir=work.dir)

    modules.gene.centers1 <- do.call(rbind, lapply(modules.gene.centers, function(x) { head(x, n=1) }))

    if (TRUE) {
        png(file=sprintf("%s/centroids.gene.final.png", work.dir), width=1000, height=500)
        pheatmap(
            modules.gene.centers1,
            cluster_rows=T, cluster_cols=T,
            show_colnames=show_colnames)
        dev.off()
    }

    res
}

makeNetworkWithCorTest <- function(gene.cor.test, base, network, eps=0) {
    gene.cor.pval <- sapply(gene.cor.test, function(x) x$p.value)
    gene.cor.rho <- sapply(gene.cor.test, function(x) unname(x$estimate))
    fake.gene.de <- data.frame(
        ID=names(gene.cor.pval),
        pval=gene.cor.pval,
        log2FC=gene.cor.rho)
    fake.gene.de$logPval <- log2(fake.gene.de$pval)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    es.rn <- makeExperimentSet(network, gene.de=fake.gene.de, reactions.as.edges=F, plot=F)
    rxn.fdr <- 10^-base
    es.rn.scored <- scoreNetwork(es.rn, rxn.fdr=rxn.fdr)
    subnet.scored <- es.rn.scored$subnet.scored

    #V(subnet.scored)$score <- -0.1
    #V(subnet.scored)$score <- -0.01
    V(subnet.scored)[nodeType == "met"]$score <- eps

    repeats <- table(V(subnet.scored)[score > 0]$origin)
    repeats <- names(repeats[repeats >= 4])
    V(subnet.scored)[origin %in% repeats]$score <- eps

    es.rn.scored$subnet.scored <- subnet.scored
    es.rn.scored
}

makeAtomNetworkWithCorTest <- function(gene.cor.test, base, eps=0) {
    gene.cor.pval <- sapply(gene.cor.test, function(x) x$p.value)
    gene.cor.rho <- sapply(gene.cor.test, function(x) unname(x$estimate))
    fake.gene.de <- data.frame(
        ID=names(gene.cor.pval),
        pval=gene.cor.pval,
        log2FC=gene.cor.rho)
    fake.gene.de$logPval <- log2(fake.gene.de$pval)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    net <- makeNetwork(gene.de = fake.gene.de)
    fb <- fitBumModel(gene.cor.pval, plot = F)
    rxn.fdr <- 10^-base
    gene.scores <- GAM:::scoreValue(fb, gene.cor.pval, rxn.fdr)
    V(net)$score <- eps
    E(net)$score <- gene.scores[E(net)$origin]
    net
}

makeENetworkWithCorTest <- function(gene.cor, gene.cor.pval, base, network, eps=0) {
#     gene.cor.pval <- sapply(gene.cor.test, function(x) x$p.value)
#     gene.cor.rho <- sapply(gene.cor.test, function(x) unname(x$estimate))
    fake.gene.de <- data.frame(
        ID=names(gene.cor),
        pval=gene.cor.pval,
        log2FC=gene.cor)
    fake.gene.de$logPval <- log2(fake.gene.de$pval)
    fake.gene.de <- merge(fake.gene.de,
                          unique(reflink[, list(Entrez, symbol)]),
                          by.x="ID", by.y="Entrez", all.x=T)

    fake.gene.de <- fake.gene.de[fake.gene.de$log2FC > 0, ]

    es.re <- makeExperimentSet(network, gene.de=fake.gene.de, reactions.as.edges=T, plot=F)
    net <- es.re$subnet
    rxn.fdr <- 10^-base
    gene.cor.fb <- fitBumModel(gene.cor.pval[gene.cor > 0], plot=F)

    if (gene.cor.fb$lambda > 1-1e-3) {
        logwarn("lambda == 1")
        return(NULL)
    }
    gene.scores <- GAM:::scoreValue(gene.cor.fb, gene.cor.pval[gene.cor > 0], rxn.fdr)

    V(net)$score <- eps
    E(net)$score <- gene.scores[E(net)$origin]
    # E(net)[log2FC < 0]$score <- -10000
    es.re$subnet.scored <- net
    es.re
}

gamCluster2 <- function(gene.exprs, tag, work.dir, network=kegg.mouse.network, show_colnames=T, cor.method="pearson") {
    stopifnot(require(cluster))
    stopifnot(require(GAM))
    stopifnot(require(igraph))
    stopifnot(require(pheatmap))
    stopifnot(require(data.table))
    stopifnot(require(BioNet))
    dir.create(work.dir, showWarnings=F, recursive=T)
    metGenes <- unique(network$rxn2gene$gene)
    gene.exprs <- gene.exprs[rownames(gene.exprs) %in% metGenes, ]
    #gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$gene.id.map$Entrez, ]
    gene.exprs <- gene.exprs[apply(gene.exprs, 1, sd) > 0,]

    #gene.pairwise.cor <- cor(t(gene.exprs))
    gene.pairwise.cor <- cor(t(gene.exprs), method=cor.method)
    diag(gene.pairwise.cor) <- 0
    gene.cor.dist <- as.dist(1 - gene.pairwise.cor)
    gene.pairwise.cor.sd <- sd(gene.pairwise.cor[upper.tri(gene.pairwise.cor)])

    gK <- 2
    gK <- 10
    reorder <- sample(gK)
    set.seed(42)
    gene.pam <- pam(gene.cor.dist, k=gK)
    clustersOrder <- order(table(gene.pam$clustering ), decreasing=T)
    gene.pam$clustering <- match(gene.pam$clustering, clustersOrder)
    names(gene.pam$clustering) <- rownames(gene.exprs)

    gene.pam$medoids <- gene.pam$medoids[clustersOrder]
    gene.pam$clusinfo <- gene.pam$clusinfo[clustersOrder, ]
    if (FALSE) {
        pheatmap(
            normalize.rows(gene.exprs[order(reorder[gene.pam$clustering]), ]),
            cluster_rows=F, cluster_cols=T,
            show_rownames=F, show_colnames=show_colnames)
        pheatmap(
            normalize.rows(gene.exprs[gene.pam$medoids, ]),
            cluster_rows=F, cluster_cols=T,
            show_rownames=T, show_colnames=show_colnames)

        pheatmap(normalize.rows(gene.exprs),
            cluster_cols = T,
            cluster_rows = T,kmeans_k = gK,
            show_rownames=T, show_colnames=show_colnames)
        tt <- gene.exprs[gene.pam$medoids, ]
        rownames(tt) <- paste0(rownames(tt), "_", gene.pam$clusinfo[, "size"])
        pheatmap(
            normalize.rows(rbind(cc, tt)),
            cluster_rows=F, cluster_cols=T,
            show_rownames=T, show_colnames=show_colnames)

    }

    # gene.dist.to.clusters <- as.matrix(gene.cor.dist)[, gene.pam$medoids]
    # gene.dist.to.cluster.closest <- apply(gene.dist.to.clusters, 1, min)
    # gene.dist.to.cluster.which <- apply(gene.dist.to.clusters, 1, which.min)

    if (FALSE) {
        z <- sample(gene.pairwise.cor[upper.tri(gene.pairwise.cor)], 10000)
        qplot(z, binwidth=1/64)
        qplot(sample=z, stat="qq", dparams=list(sd=sd(z))) + coord_equal(ratio=1) + xlim(-1, 1) + ylim(-1, 1)
    }

    modules.gene <- list()
    modules.gene.centers <- list()

    gene.bases <- c(rep(0, gK))
    baseMean <- apply(gene.exprs, 1, mean)

    i <- 1
    
    #dir.create(sprintf("%s/refined", work.dir), showWarnings=FALSE, recursive=TRUE)
    
    idxs <- 6:10
    idxs <- c(2, 7)
    idxs <- 4
    #idxs <- c(3, 1:2, 4:10)
    idxs <- seq_len(gK)

    for (i in idxs) {
    #results <- mclapply(idx, mc.cores = 1, FUN = function(i) {
        stepBig <- 1
        step <- 0.25
        cluster.genes <- names(gene.pam$clustering)[gene.pam$clustering == i]
        cluster.genes.neg <- NULL
        #weights <- rep(1, length(cluster.genes))

        medoid <- gene.exprs[gene.pam$medoids[i],]
        medoid <- (medoid - mean(medoid)) / sd(medoid)

        it = 0
        base <- gene.bases[i]


        centers <- rbind(medoid)
        rownames(centers)[1] <- paste0("medoid.", i)

        redo <- TRUE
        while (TRUE) {
            center <- getCenter(gene.exprs, cluster.genes, cluster.genes.neg, method=cor.method)

            cor.diff <- 1 - max(apply(centers, 1, cor, y=center, method=cor.method))

            centers <- rbind(center, centers)
            rownames(centers)[1] <- paste0("center.", i, ".r", it)
            if (T) {
                pheatmap(
                    centers[, order(center)],
                    cluster_rows=F, cluster_cols=F,
                    show_colnames=show_colnames)
            }

            if (!redo && cor.diff < 0.001) {
                break
            }
            redo <- FALSE

            loginfo("Difference: %s", cor.diff)
            loginfo("Core: %s", paste0(cluster.genes, collapse=" "))
            loginfo("iteration #%s, base=%s", it, base)

            gene.cor <- cor(t(gene.exprs), as.matrix(center), method=cor.method)[,1]

#             gene.cor.test <-
#                     apply(gene.exprs, 1, cor.test, y=center, exact=F, alternative="t", method="spearman")
#             gene.cor.pval <- sapply(gene.cor.test, function(x) x$p.value)
#             gene.cor.rho <- sapply(gene.cor.test, function(x) unname(x$estimate))

            gene.cor.pval <- 2*pnorm(-abs(gene.cor), mean = 0, sd=gene.pairwise.cor.sd)
            if (base <= 0) {
                gene.cor.fb <- fitBumModel(gene.cor.pval, plot = T)

                if (gene.cor.fb$lambda > 1-1e-3) {
                    module <- NULL
                    logwarn("Lambda == 1 for cluster %s", i)
                    break
                }

                base <- step * ceiling(
                    -log10(GAM:::recommendedFDR(gene.cor.fb, gene.cor.pval, num.positive=30)) / step
                )

                if (base <= 1) {
                    module <- NULL
                    logwarn("Base <= 1 for cluster %s", i)
                    break
                }

                loginfo("iteration #%s, base=%s", it, base)
            }

            # es.rn.scored <- makeNetworkWithCorTest(gene.cor.test, base, network)
            es.re.scored <- makeENetworkWithCorTest(gene.cor, gene.cor.pval, base, network)

            module <- findModule(es.re.scored, solver.gmwcs.big)

            module

            while (is.null(module)) {
                Sys.sleep(30)
                #module <- findModule(es.re.scored, solver.gmwcs)
                module <- findModule(es.re.scored, solver.gmwcs.big)
            }
            loginfo("Solved to optimality: %s", attr(module, "optimal"))

            if (F) {
                dir.create(sprintf("%s/refined/", work.dir), showWarnings = F)
                saveModuleToXgmml(
                    module,
                    file=sprintf("%s/refined/gene.%s.r%s.xgmml", work.dir, i, it),
                    name=sprintf("%s.gene.%s.refined.r%s.b%s", tag, i, it, base))

            }

            modules.gene[[i]] <- module
            modules.gene.centers[[i]] <- centers

            #V(module)$betweenness <- betweenness(module)
            #bet.median <- median(V(module)$betweenness)
            #core <- unique(na.omit(V(module)[betweenness < bet.median]$origin))
            core <- unique(na.omit(E(module)[score > 0 & log2FC > 0]$origin))
            core.neg <- unique(na.omit(E(module)[score > 0 & log2FC < 0]$origin))
            core.neg <- c()
            core.size <- length(core) + length(core.neg)

            # ???
            cluster.genes.new <- intersect(core, rownames(gene.exprs))
            cluster.genes.new.neg <- intersect(core.neg, rownames(gene.exprs))
            if (attr(module, "optimal") && core.size <= 100 && length(cluster.genes.new) + length(cluster.genes.new.neg) >= 10) {
                cluster.genes <- cluster.genes.new
                cluster.genes.neg <- cluster.genes.new.neg
            }

            #         cluster.genes<- intersect(E(module)$RefSeq, rownames(gene.exprs))
            #             saveModuleToXgmml(
            #                 module,
            #                 file=sprintf("%s/refined/c%s.gene.%s.r%s.xgmml", work.dir, gK, i, it),
            #                 name=sprintf("%s.gene.%s.refined.r%s.b%s", tag, i, it, base))

            if (core.size > 100 || !attr(module, "optimal")) {
                base <- base + stepBig
                redo <- TRUE
                step <- step / 2
                loginfo("Big module, decreasing step to %s", step)
            } else if (core.size < 30 || get.score(module) <= 0) {
                base <- base - step
                redo <- TRUE
            }

            if (base <= 0) {
                module <- NULL
                logwarn("Base <= 0 for cluster %s", i)
                break
            }
            it = it + 1
        }

        if (!is.null(module)) {
            es.re.scored <- makeENetworkWithCorTest(gene.cor, gene.cor.pval, base, network, eps=-0.01)
            #         stopifnot(sum(V(es.rn.scored$subnet.scored)$score) < 0)
            #
            #
            #
            #         #module <- findModule(es.rn.scored, solver.mc)
            module <- findModule(es.re.scored, solver.gmwcs.big)
            while (is.null(module)) {
                Sys.sleep(30)
                module <- findModule(es.re.scored, solver.gmwcs.big)
                #module <- findModule(es.re.scored, solver.gmwcs)
            }

            modules.gene[[i]] <- module

            #

            #         pheatmap(
            #             centers,
            #             cluster_rows=F, cluster_cols=T,
            #             show_colnames=F)
            #
            saveModuleToXgmml(
                module,
                file=sprintf("%s/c%s.gene.%s.xgmml", work.dir, gK, i),
                name=sprintf("%s.gene.%s.refined.b%s", tag, i, base))
        }
        #list(module=module, center=center)
    }

    #     pheatmap(
    #         do.call(rbind, modules.gene.centers),
    #         cluster_rows=F, cluster_cols=T,
    #         file=sprintf("%s/centroids.gene.all.png", work.dir))
    #
    res <- list(modules=modules.gene, centers=modules.gene.centers)
    writeClusterGenesE(exprs=gene.exprs, cluster=res, tag=tag, work.dir=work.dir)

    modules.gene.centers1 <- do.call(rbind, lapply(modules.gene.centers, function(x) { head(x, n=1) }))

    if (TRUE) {
        png(file=sprintf("%s/centroids.gene.final.png", work.dir), width=1000, height=500)
        pheatmap(
            modules.gene.centers1,
            cluster_rows=T, cluster_cols=T,
            show_colnames=show_colnames)
        dev.off()
    }

    res
}

writeClusterGenes <- function(exprs, cluster, tag, work.dir) {
    for (i in seq_along(cluster$centers)) {
        if (is.null(dim(cluster$centers[[i]]))) {
            ci <- cluster$centers[[i]]
        } else {
            ci <- cluster$centers[[i]][1,]
        }

        ci.cor <- apply(exprs, 1, cor, y=ci)
        ci.cor <- ci.cor[which(ci.cor > 0.5)]
        ci.table <- data.table(gene=names(ci.cor), cor=ci.cor, symbol=reflink[match(names(ci.cor), Entrez), symbol])
        ci.table <- ci.table[order(cor, decreasing=TRUE), ]
        write.tsv(ci.table, file=file.path(work.dir, sprintf("%s.c%s.tsv", tag, i)))
        m <- cluster$modules[[i]]
        writeLines(unique(V(m)[nodeType == "rxn"]$symbol), file.path(work.dir, sprintf("%s.c%s.module.genes.tsv", tag, i)))
        writeLines(unique(V(m)[nodeType == "met"]$label), file.path(work.dir, sprintf("%s.c%s.module.mets.tsv", tag, i)))
    }
}

writeClusterGenesE <- function(exprs, cluster, tag, work.dir, cor.method="pearson") {
    gene.pairwise.cor <- cor(t(exprs), method=cor.method)
    diag(gene.pairwise.cor) <- 0
    gene.pairwise.cor.sd <- sd(gene.pairwise.cor[upper.tri(gene.pairwise.cor)])



    for (i in seq_along(cluster$centers)) {
        if (is.null(cluster$centers[[i]])) {
            next
        }

        if (is.null(dim(cluster$centers[[i]]))) {
            ci <- cluster$centers[[i]]
        } else {
            ci <- cluster$centers[[i]][1,]
        }

        ci.cor <- apply(exprs, 1, cor, y=ci, method=cor.method)
        ci.cor <- ci.cor[which(ci.cor > 0.5)]
        cor.pval <- 2*pnorm(-abs(ci.cor), mean = 0, sd=gene.pairwise.cor.sd)


        ci.table <- data.table(gene=names(ci.cor), cor=ci.cor, symbol=reflink[match(names(ci.cor), Entrez), symbol], pval=cor.pval)
        ci.table <- ci.table[order(cor, decreasing=TRUE), ]
        write.tsv(ci.table, file=file.path(work.dir, sprintf("%s.c%s.tsv", tag, i)))
        m <- cluster$modules[[i]]
        writeLines(unique(E(m)$symbol), file.path(work.dir, sprintf("%s.c%s.module.genes.tsv", tag, i)))
        writeLines(unique(E(m)$label), file.path(work.dir, sprintf("%s.c%s.module.mets.tsv", tag, i)))
    }
}

geneSimilarity <- function(module1, module2) {
    setSimilarity(V(module1)[nodeType == "rxn"]$origin,
                  V(module2)[nodeType == "rxn"]$origin)
}

loadClusters <- function(dir) {
    files <- list.files(dir, pattern="c\\d+\\.tsv$")
    tags <- gsub(".tsv$", "", files)
    gene.files <- file.path(dir, paste0(tags, ".module.genes.tsv"))
    genes <- lapply(gene.files, function(f) na.omit(read.table(f, stringsAsFactors=F)[[1]]))
    met.files <- file.path(dir, paste0(tags, ".module.mets.tsv"))
    mets <- lapply(met.files, function(f) na.omit(read.table(f, stringsAsFactors=F, quote="", sep="\t")[[1]]))
    clusters <- list()
    for (i in seq_along(tags)) {
        tag <- tags[[i]]
        clusters[[tag]] <- list(genes=genes[[i]], mets=mets[[i]])
    }
    clusters
}

adjustClusters <- function(gene.exprs, clusters, tag, work.dir, network, base=0.5, cut=base, show_colnames=T) {
    dir.create(work.dir, recursive=T, showWarnings=F)
    gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$gene.id.map$Entrez, ]
    gene.exprs <- gene.exprs[apply(gene.exprs, 1, sd) > 0,]

    mergeModules <- function(genes, groups) {
        centers <- list()
        modules <- list()
        for (i in seq_along(groups)) {
            message(sprintf("group #%s", i))
            cluster.genes <- unlist(genes[groups[[i]]])
            center <- getCenter(gene.exprs, cluster.genes)

            centers[[i]] <- center

            gene.dist.to.cluster <- 1 - cor(t(gene.exprs), as.matrix(center))[,1]
            es.rn.scored <- makeNetworkWithDist(gene.dist.to.cluster, base, network)

            #module <- findModule(es.rn.scored, solver.mc)
            module <- findModule(es.rn.scored, heinz2.solver("/usr/local/lib/heinz2/heinz", nthreads=6, timeLimit=30))

            modules[[i]] <- module
            module
        }
        list(modules=modules, centers=centers)
    }

    if ("modules" %in% names(clusters)) {
        genes <- lapply(clusters$modules, function(m) unique(V(m)[nodeType=="rxn"]$origin))
    } else {
        genes <- lapply(clusters, function(c) network$gene.id.map[na.omit(match(unique(c$genes), Symbol)), Entrez] )
    }

    res <- mergeModules(genes, groups=as.list(seq_along(genes)))

    while (TRUE) {
        centers <- lapply(genes, getCenter, gene.exprs=gene.exprs)
        h <- hclust(as.dist(1-pairwiseCompare(cor, centers)), method="average")
        if (FALSE) {
            plot(h)
        }

        groups <- cutree(h, h=cut)
        groups <- split(seq_along(groups), groups)

#         mm <- pairwiseCompare(setSimilarity, genes)
#         mm <- pairwiseCompare(cor, centers)
#         if (FALSE) {
#             pheatmap(mm)
#         }
#         groups <- igraph::clusters( graph.adjacency(mm >= 0.5))

#         if (max(groups$csize) == 1 & !is.null(res)) {
#             break
#         }
#         groups <- split(seq_along(groups$membership), groups$membership)

        groups <- groups[sapply(groups, length) > 1]
        if (length(groups) == 0) {
            break
        }
        res1 <- mergeModules(genes, groups)
        res$modules[unlist(groups)] <- NULL
        res$centers[unlist(groups)] <- NULL
        res$modules <- c(res$modules, res1$modules)
        res$centers <- c(res$centers, res1$centers)

        genes <- lapply(res$modules, function(m) unique(V(m)[nodeType=="rxn"]$origin))
    }

    toKeep <- sapply(res$modules, function(m) sum(V(m)$score > 0) >= 5)
    res$modules <- res$modules[toKeep]
    res$centers <- res$centers[toKeep]


    for (i in seq_along(res$modules)) {
        saveModuleToXgmml(
            res$modules[[i]],
            file=sprintf("%s/c%s.xgmml", work.dir, i),
            name=sprintf("%s.gene.%s.refined", tag, i))
    }

    png(file=sprintf("%s/clusters.final.png", work.dir), width=1000, height=500)

    pheatmap(
        normalize.rows(gene.exprs[unlist(genes), ]),
        cluster_rows=F, cluster_cols=T,
        show_colnames=show_colnames)
    dev.off()

    writeClusterGenes(exprs=gene.exprs, cluster=res, tag=tag, work.dir=work.dir)
    res
}

adjustClusters2 <- function(gene.exprs, clusters, tag, work.dir, network, base, cut=0.5, show_colnames=T, cor.method="pearson") {
    dir.create(work.dir, recursive=T, showWarnings=F)
    gene.exprs <- gene.exprs[rownames(gene.exprs) %in% network$rxn2gene$gene, ]
    gene.exprs <- gene.exprs[apply(gene.exprs, 1, sd) > 0,]

    gene.pairwise.cor <- cor(t(gene.exprs), method=cor.method)
    diag(gene.pairwise.cor) <- 0
    gene.pairwise.cor.sd <- sd(gene.pairwise.cor[upper.tri(gene.pairwise.cor)])

    mergeModules <- function(genes, groups) {
        centers <- list()
        modules <- list()
        for (i in seq_along(groups)) {
            loginfo("group #%s", i)
            cluster.genes <- unlist(genes[groups[[i]]])
            center <- getCenter(gene.exprs, cluster.genes, method = cor.method)

            centers[[i]] <- center

            gene.cor <- cor(t(gene.exprs), as.matrix(center), method=cor.method)[,1]
            gene.cor.pval <- 2*pnorm(-abs(gene.cor), mean = 0, sd=gene.pairwise.cor.sd)

            es.re.scored <- makeENetworkWithCorTest(gene.cor, gene.cor.pval,
                                                    base, network, eps=-0.01)

            if (is.null(es.re.scored)) {
                modules[[i]] <- NULL
                logwarn("No network for core %s", paste0(cluster.genes, collaspe=" "))
                next
            }

            module <- findModule(es.re.scored, solver.gmwcs.big)
            while (is.null(module)) {
                Sys.sleep(30)
                module <- findModule(es.re.scored, solver.gmwcs.big)
                #module <- findModule(es.re.scored, solver.gmwcs)
            }

            loginfo("Solved to optimality: %s", attr(module, "optimal"))

            modules[[i]] <- module
            module
        }

        res <- list(modules=modules, centers=centers)
        res$centers[sapply(res$modules, is.null)] <- NULL
        res$modules[sapply(res$modules, is.null)] <- NULL
        res
    }

    clusters$centers[sapply(clusters$modules, is.null)] <- NULL
    clusters$modules[sapply(clusters$modules, is.null)] <- NULL

    genes <- lapply(clusters$modules, function(m) unique(E(m)[log2FC > 0 & score > 0]$origin))

    res <- mergeModules(genes, groups=as.list(seq_along(genes)))

    while (TRUE) {
        centers <- lapply(genes, getCenter, gene.exprs=gene.exprs)
        h <- hclust(as.dist(1-pairwiseCompare(cor, centers)), method="average")
        if (FALSE) {
            plot(h)
        }

        groups <- cutree(h, h=cut)
        groups <- split(seq_along(groups), groups)

        #         mm <- pairwiseCompare(setSimilarity, genes)
        #         mm <- pairwiseCompare(cor, centers)
        #         if (FALSE) {
        #             pheatmap(mm)
        #         }
        #         groups <- igraph::clusters( graph.adjacency(mm >= 0.5))

        #         if (max(groups$csize) == 1 & !is.null(res)) {
        #             break
        #         }
        #         groups <- split(seq_along(groups$membership), groups$membership)

        groups <- groups[sapply(groups, length) > 1]
        if (length(groups) == 0) {
            break
        }
        res1 <- mergeModules(genes, groups)
        res$modules[unlist(groups)] <- NULL
        res$centers[unlist(groups)] <- NULL
        res$modules <- c(res$modules, res1$modules)
        res$centers <- c(res$centers, res1$centers)

        genes <- lapply(res$modules, function(m) unique(E(m)[log2FC > 0 & score > 0]$origin))
    }

    toKeep <- sapply(res$modules, function(m) length(unique(E(m)[score > 0]$origin)) >= 5)
    res$modules <- res$modules[toKeep]
    res$centers <- res$centers[toKeep]

    for (i in seq_along(res$modules)) {
        saveModuleToXgmml(
            res$modules[[i]],
            file=sprintf("%s/c%02d.xgmml", work.dir, i),
            name=sprintf("%s.gene.c%02d", tag, i))
        t <- res$modules[[i]]
        E(t)$log2FC <- E(t)$log2FC * 2
        saveModuleToDot(t, file=sprintf("%s/c%02d.dot", work.dir, i))
        rm(t)
        system2("neato", c("-Tpdf", "-O", sprintf("%s/c%02d.dot", work.dir, i)), stderr=NULL)
    }

    annotation_row <- do.call(rbind, lapply(seq_along(genes), function(i) {
        data.frame(gene=genes[[i]], cluster=sprintf("c%02d", i))
    }))

    annotation_row <- aggregate(cluster ~ gene, annotation_row, paste0, collapse="&")

    rownames(annotation_row) <- annotation_row$gene
    annotation_row$gene <- NULL

    pheatmap(
        normalize.rows(gene.exprs[unlist(genes), ]),
        file=sprintf("%s/clusters.final.png", work.dir), width=12, height=8,
        cluster_rows=F, cluster_cols=T,
        annotation_row=annotation_row,
        show_colnames=show_colnames)

    if (F){
        pheatmap(
            normalize.rows(gene.exprs[unlist(genes), ]),
            cluster_rows=F, cluster_cols=T,
            annotation_row=annotation_row,
            show_colnames=show_colnames)
    }


    writeClusterGenesE(exprs=gene.exprs, cluster=res, tag=tag,
                       work.dir=work.dir, cor.method=cor.method)
    res
}

#     clustDist <- function(genes1, genes2) {
#         center1 <- getCenter(gene.exprs, genes1)
#         center2 <- getCenter(gene.exprs, genes2)
#         1 - cor(center1, center2)
#     }
#     myHClust <- function(genes) {
#         genes1 <- genes
#         centers1 <- lapply(genes, getCenter, gene.exprs=gene.exprs)
#         for (t in seq_len(length(genes) - 1)) {
#             mm <- pairwiseCompare(cor, centers1)
#             pheatmap(mm)
#             diag(mm) <- NA
#             k <- arrayInd(which.max(mm), .dim=dim(mm))
#             i <- k[,1]
#             j <- k[,2]
#             genes1[[i]] <- c(genes1[[i]], genes1[[j]])
#             genes1[[j]] <- NULL
#             centers1[[i]] <- getCenter(gene.exprs, genes1[[i]])
#             centers1[[j]] <- NULL
#         }
#     }

gmwcs.solver2 <- function (gmwcs, nthreads = 1, timeLimit = -1, nodes.group.by=NULL, edges.group.by=NULL, c.size=10000, group.only.positive=F, quiet=T) {
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        synonyms.file <- file.path(graph.dir, "synonyms.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (!score.edges) {
            E(network)$score <- 0
        }

        synonyms <- c()

        nt <- get.vertex.attributes(network)
        if (!is.null(nodes.group.by)) {
            synonyms <- c(synonyms, aggregate(as.formula(sprintf("name ~ %s", nodes.group.by)),
                                              data=nt, paste0, collapse=" ")$name)
        }
        nt <- rename(nt[, c("name", "score")], c("name"="#name"))

        et <- get.edge.attributes(network, include.ends = T)
        if (!is.null(edges.group.by)) {
            etx <- if (group.only.positive) { et[et$score > 0,] } else { et }
            if (nrow(etx) > 0) {
                synonyms <- c(synonyms,
                              aggregate(name ~ edges.group.by,
                                        data=list(
                                            name=sprintf("%s -- %s", etx$from, etx$to),
                                            edges.group.by=etx[[edges.group.by]]),
                                        paste0, collapse=" ")$name
                )
            }
        }
        et <- rename(et[, c("from", "to", "score")], c("from"="#from"))

        write.tsv(nt, file=nodes.file)
        write.tsv(et, file=edges.file)
        writeLines(sprintf("%s\n", synonyms), con=synonyms.file)


        out <- system2(gmwcs, c("--nodes", nodes.file,
                         "--edges", edges.file,
                         if (length(synonyms) > 0) c("--synonyms", synonyms.file, "-B", 1) else NULL,
                         "--threads", nthreads,
                         "--timelimit", timeLimit,
                         "-c", c.size,
                         "--break"
        ), stdout=if (quiet) T else "")
        solution.file <- paste0(nodes.file, ".out")
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- GAM:::readGraph(node.file = solution.file,
                               edge.file = paste0(edges.file, ".out"),
                               network = network)
        if (quiet) {
            attr(res, "optimal") <- any(grepl("SOLVED TO OPTIMALITY", out))
        }
        return(res)
    }
}

sgmwcs.batchSolver <- function(gmwcs, nthreads = 1, timeLimit = -1, nodes.group.by=NULL, edges.group.by=NULL, c.size=NULL, group.only.positive=F, quite
                               =T) {
    function(networks) {
        N <- length(networks)
        graphBatch.dir <- tempfile("graphBatch")
        dir.create(graphBatch.dir)
        graph.dirs <- file.path(graphBatch.dir, seq_len(N))

        instances <- lapply(seq_len(N), function(i) {
            network <- networks[[i]]
            score.edges <- "score" %in% list.edge.attributes(network)
            score.nodes <- "score" %in% list.vertex.attributes(network)
            if (!score.nodes) {
                V(network)$score <- 0
            }
            if (!score.edges) {
                E(network)$score <- 0
            }

            graph.dir <- graph.dirs[i]
            dir.create(graph.dir)

            writeGmwcsInstance(graph.dir = graph.dir,
                               network=network,
                               nodes.group.by = nodes.group.by,
                               edges.group.by = edges.group.by,
                               group.only.positive = group.only.positive
            )
        })

        system2(gmwcs, c("-a", sprintf("1-%s", N),
                         "-d", graphBatch.dir,
                         "--threads", nthreads,
                         "--timelimit", timeLimit,
                         if (!is.null(c.size)) c("-c", c.size) else NULL,
                         "--break"
        ))

        res <- lapply(seq_len(N), function(i) {
            network <- networks[[i]]
            instance <- instances[[i]]
            graph.dir <- graph.dirs[i]
            solution.file <- paste0(instance$nodes.file, ".out")

            if (!file.exists(solution.file)) {
                warningf("Solution file '%s' not found", solution.file)
                NULL
            } else {
                GAM:::readGraph(node.file = solution.file,
                                edge.file = paste0(instance$edges.file, ".out"),
                                network = network)
            }
        })
    }
}

dualGraph <- function(m, what) {
    t1 <- get.edge.attributes(m, include.ends = T)[, c("from", "to", what)]
    t2 <- rename(t1, c("to"="from", "from"="to"))
    tt <- unique(rbind(t1, t2)[, c("from", what)])
    dm <- graph.from.tables(edge.table=unique(merge(tt, tt, by="from")
                                              [, c(paste0(what, ".x"),
                                                   paste0(what, ".y"))]))
    dm
}

makeENetworkWithScore <- function(score, base, network) {
    fake.gene.de <- data.frame(
        ID=names(score),
        pval=2 ^ -(score),
        score=score,
        log2FC=score)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    suppressWarnings(es.re <- makeExperimentSet(network, gene.de=fake.gene.de, reactions.as.edges=T, plot=F))
    subnet.scored <- es.re$subnet

    #V(subnet.scored)$score <- -0.1
    V(subnet.scored)$score <- 0
    rxn.vertices <- V(subnet.scored)[nodeType == "rxn"]
    #     E(subnet.scored)$score <-
    #         -(1 - pmax(1 - E(subnet.scored)$dist.closest, 0)) +
    #         base
    #

    es.re.scored <- es.re
    es.re.scored$subnet.scored <- subnet.scored
    es.re.scored
}

makeANetworkWithScore <- function(score, base, network) {
    fake.gene.de <- data.frame(
        ID=names(score),
        pval=2 ^ -(score),
        score=score,
        log2FC=score)
    fake.gene.de <- merge(fake.gene.de, unique(reflink[, list(Entrez, symbol)]), by.x="ID", by.y="Entrez", all.x=T)

    fake.gene.de$logPval <- log2(fake.gene.de$pval)


    net <- makeNetwork(gene.de = fake.gene.de)
    V(net)$score <- 0

    list(subnet.scored=net)
}