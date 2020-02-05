scNPF R package
====================

An integrative framework assisted by network propagation and network fusion for pre-processing of single-cell RNA-seq data

About
====================
scNPF is an R package for pre-processing of single-cell RNA-seq data by leveraging the context-specific topology inherent in the given data and the network information from priori gene-gene interaction networks. The basic scNPF framework consists of two modules (Figure 1), including scNPF-propagation for imputing dropouts and scNPF-fusion for fusing multiple smoothed expression matrices to learn a cell-cell similarity matrix. The scNPF framework is highly integrative and flexible in that the two modules are independent but interconnected. scNPF-propagation involves a network propagation process based on random walk with restart (RWR) on a given gene-gene interaction network to obtain a distribution for each node (gene), which captures its relevance to all other genes in the network. This process takes the global connectivity patterns of the interaction network into account for profiling the topological context of each gene. More importantly, this module contains two modes of propagation, including the priori mode that uses a publicly available interaction network and the context mode that is solely based on the given scRNA-seq data set. The output of scNPF-propagation is a propagated gene-cell expression matrix, which could be used as input for scNPF-fusion. scNPF-fusion constructs a sample-similarity network for each propagated expression matrix and then integrates different networks into a single cell-cell similarity network based on a nonlinear combination method. The learned similarity matrix from scNPF-fusion or the smoothed expression matrix from scNPF-propagation can be used as inputs for other existing scRNA-seq pipelines or tools for downstream analyses, such as cell type clustering, dimension reduction, and visualization.

![image](https://github.com/BMILAB/scNPF/blob/master/images/Schematic%20diagram%20of%20the%20scNPF%20framework.png)
                       Figure 1. Schematic diagram of the scNPF framework.

Installing scNPF
=============
Mandatory 
---------

* R (>3.1). [R 3.5](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [igraph](https://cran.r-project.org/web/packages/igraph/index.html), [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html), [foreach](https://cran.r-project.org/web/packages/foreach/index.html), [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [plyr](https://cran.r-project.org/web/packages/plyr/index.html),  

Suggested R Packages
---------
* [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html),  

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/scNPF")
library(scNPF)
```

Preparations
====================

Gene expression count matrix
---------
The input to scNPF is matrix of gene expression count. The rows correspond to genes and the columns correspond to cells. In this study, we will use the human embryonic stem cells data from [yan](http://dx.doi.org/10.1038/nsmb.2660) as example.
```
##Loading gene expression count matrix
load(system.file("data","yan.Rdata",package = "scNPF")
exp.data <- yan$data
dim(exp.data)
```

Gene-gene interaction network
---------
A gene-gene interaction network (a adjacency matrix) is used to smooth expression values in the scNPF-propagation model. If users use priori mode, they should provide gene co-expression network from publicly available database or a specific gene-gene network established on your own method. In this package, we provided three human gene-gene interaction networks from different databases, including [String](https://doi.org/10.1093/nar/gks1094), [HumanNet](http://www.functionalnet.org/humannet/about.html) and [an integrated network](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741255/). If users use context mode, a gene co-expression network is automatically generated in the scNPF.pro function.


Using scNPF
=============
In order to facilitate user understanding, we use the provided example [yan data](http://dx.doi.org/10.1038/nsmb.2660) to illustrate the standard analysis work-flow of scNPF. Please refer to the [User Guide](https://github.com/BMILAB/scNPF/tree/master/doc) for full details.

Section 1 scNPF-propagation
---------
scNPF-propagation involves a network propagation process based on RWR on a given gene-gene interaction network to obtain a distribution for each node (gene), which captures its relevance to all other genes in the network. In this step, users can use priori mode or context mode.

*Use priori mode
```
##Using String network to smooth expression values.
load(system.file("data","string.Rdata",package = "scNPF"))
string.data <- scNPF.pro(x=exp.data, network=string,nThreads=8)
```
*Use context mode
```
context.data<- scNPF.pro(x=exp.data, network="context",nThreads=8)
```
The output of function scNPF.pro is a propagated gene-cell expression matrix , which could be used as input for scNPF-fusion, and also as the input for many other single cell tools to perform downstream analyses like dimension reduction, clustering, and visualization.

Section 2 scNPF-fusion
---------
 scNPF-fusion constructs a sample-similarity network for each propagated expression matrix and then integrates these networks into a single cell-cell similarity network based on a nonlinear combination method. For example, we takes the propagated matrices from scNPF-propagation using the priori mode with the String network and the context mode as inputs and learns a matrix of similarities between cells by network fusion.
 ```
 ##Construction a cell-by-cell similarity matrix.
similarity<-scNPF.fus(data=list(string=string.data,context=context.data))
 ```

