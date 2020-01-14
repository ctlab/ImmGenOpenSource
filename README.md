# ImmGen Open Source
#### scripts for RNAseq data processing by GAM-clustering

> The Open Source [ImmGen](http://www.immgen.org/) project ([GSE122108](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122108)) is a collaborative effort devoted to RNAseq profiling of *ex vivo* sorted mononuclear phagocytes. 

> GAM-clustering provides metabolic variability within dataset using a novel network-based computational approach that utilizes cellular transcriptional profiles as proxies. The metabolic network of reactions from KEGG database is presented as a graph that has vertices corresponding to metabolites and the edges corresponding to the reactions with the expressed genes. In the graph the method tries to find a set of connected subgraphs, with each corresponding well to a certain gene expression pattern. Curret analysis reveals the major metabolic features associated with different subpopulations and highlights a number of metabolic modules that are specific to individual cell types, tissues of residence, or developmental stages.

To explore data visit the following links:
- [Gene expression heatmap](https://artyomovlab.wustl.edu/phantasus/?preloaded=ImmGen_total_Eduw0mei4)
- [PCA with samples annotation](http://artyomovlab.wustl.edu/publications/supp_materials/Immgen/PCADatasetOverview.html)
- [Heatmap of metabolic modules](https://artyomovlab.wustl.edu/phantasus/?session=x039baa087a35e7) 

## Requirments
- [R](https://www.r-project.org/)
- [sgmwcs-solver](https://github.com/ctlab/sgmwcs-solver) 	
- [KEGG mouse metabolic network](GAM) 

## Input data
[Raw counts](Data/OSMNP_unnormalized_genes_count_10_3_18.count_table) are processed by [rawDataProcessing.R]() script and the output object `es.top12k` has the following structure:

``` r
> load("Data/337_es.top12k.Rda")
> dplyr::glimpse(exprs(es.top12k))
## num [1:12000, 1:337] 14.7 12.7 13 12.1 13.4 ...
## - attr(*, "dimnames")=List of 2
##  ..$ : chr [1:12000] "Actb" "Cst3" "Fth1" "Eef1a1" ...
##  ..$ : chr [1:337] "MF.64pLYVEpIIn.Ao.1" "MF.64pLYVEpIIn.Ao.2" "MF.64pLYVEpIIn.Ao.3" ...
> Biobase::exprs(es.top12k)[1:3,1:3]
##        MF.64pLYVEpIIn.Ao.1 MF.64pLYVEpIIn.Ao.2 MF.64pLYVEpIIn.Ao.3
## Actb              14.71612            15.00513            14.76655
## Cst3              12.67268            12.83560            12.58577
## Fth1              12.95467            13.27726            13.06159
```

## Modules deriving
The initial patterns are defined using k-means clustering on gene expression matrix and then are refined in an iterative process using the network connections ([modulesDeriving.R]()).
The final output presents a set of specific subnetworks (also called metabolic modules) that reflect metabolic variability within a given transcriptional dataset. 
Each metabolic module is a piece of metabolic network whose gene expression has correlated expression pattern across all dataset. The following graph and heatmap represent network and constituting genes' expression for module 5, correspondingly:
![module5](/readmePics/github.pic.m5.png "network and gene expression heatmap for module 5")
Averaged gene expression of all modules is represented at the following summary heatmap:
![centers](/readmePics/github.m.centers.png "centers heatmap")

## Modules annotation
Functional annotation of obtained modules is based on KEGG and Reactome canonical pathways ([modulesAnnotation.R]()).
The following example is devoted to module 5 (k - number of module genes in a particular pathway, K - number of genen in a particular pathway):
``` r
> paths <- data.table::fread("Data/m.5.pathways_mod.tsv")
> paths[1:3,]
##           PATHID         pval  k   K         padj                 PATHNAME				genes
## 1:  R-MMU-191273 1.009604e-48 17  24 1.237774e-45 Cholesterol biosynthesis	Hmgcs1 Hmgcr Msmo1 Cyp51 Mvd ...
## 2: R-MMU-8957322 9.801101e-39 17  67 6.008075e-36   Metabolism of steroids	Hmgcs1 Hmgcr Msmo1 Cyp51 Mvd ...
## 3:  R-MMU-556833 1.406903e-27 18 395 5.749543e-25     Metabolism of lipids	Hmgcs1 Hmgcr Msmo1 Cyp51 Aacs ...
```