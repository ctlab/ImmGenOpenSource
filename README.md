# ImmGenOpenSource
#### scripts for RNAseq data processing by GAM-clustering

Similar to GAM metabolic network analysis [[1](https://www.ncbi.nlm.nih.gov/pubmed/27098040)], a metabolic network of reactions from KEGG database is considered. The metabolic network is presented as a graph that has vertices corresponding to metabolites and the edges corresponding to the reactions with the expressed genes. In the graph the method tries to find a set of connected subgraphs, with each corresponding well to a certain gene expression pattern. The initial patterns are defined using k-means clustering on gene expression matrix and then are refined in an iterative process using the network connections.

Interactive heatmap of the whole dataset is available [here](https://artyomovlab.wustl.edu/phantasus/?preloaded=ImmGen_total_Eduw0mei4). 
Interactive heatmap of modules is available here [[2](https://artyomovlab.wustl.edu/phantasus/?preloaded=modules_ohmo8aeLj),[3](https://artyomovlab.wustl.edu/phantasus/?session=x039baa087a35e7)] and as a session (see SupplData/modules_ohmo8aeLj.json) for [Phantasus](https://artyomovlab.wustl.edu/phantasus/).

Before the analysis make sure there are [gmwcs](https://github.com/ctlab/gmwcs-solver)- and [sgmwcs](https://github.com/ctlab/sgmwcs-solver)-solvers installed. 	

