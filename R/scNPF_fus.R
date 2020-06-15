#' @title  Similarity matrix Construction
#'
#' @description constructs a sample-similarity network
#'      for each propagated expression matrix and then
#'      integrates these networks into a single cell-cell
#'      similarity network based on a nonlinear combination
#'      method.
#'
#' @param data A list stored network-smoothed gene expression
#'        matrix with two length.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param alpah Variance for local model.
#' @param T Number of iterations for the diffusion process.
#' @return A consensus similarity matrix.
#' @examples
#' ##Loading example data
#' load(system.file("data","yan.Rdata",package = "scNPF"))
#' exp.data <- yan$data[1:2000,]
#'
#' ##Smoothed gene expression based on priori mode
#' ##Using String network.
#' load(system.file("data","string.Rdata",package = "scNPF"))
#' string.data <- scNPF.pro(x=exp.data, network=string)
#'
#' ##Smmothed gene expression based on context mode
#' context.data<- scNPF.pro(x=exp.data, network="context",qt.gene=0.9,qt.cell=0.9)
#'
#' ##Construction a cell-by-cell similarity matrix.
#' similarity <- scNPF.fus(data=list(string=string.data,context=context.data))
#'
#' @export
scNPF.fus <- function(data,K=20,alpha = 0.5, T=10){
  if(is.null(names(data))){
    names(data)<-c(1:length(data))
  }
  data.names <- names(data)
  data.dist <- list()
  data.simi <- list()
  for(i in 1:length(data.names)){
    data.temp <- data[[data.names[i]]]
    #calculate pearson distance
    dist.temp <- as.matrix(1-WGCNA::cor(data.temp, method = "pearson"))
    data.dist[[data.names[i]]] <- dist.temp
    data.simi[[data.names[i]]] <- affinityMatrix(Diff=dist.temp, K=K, sigma=alpha)
  }
  W <- SNF(Wall=data.simi,K=K,t=T)
  return(W)
}



