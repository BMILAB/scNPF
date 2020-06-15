#' @title  Imputing dropout values of scRNA-seq data.
#'
#' @description scNPF-propagation for imputing dropouts and correcting expression
#'         expression measurements.scNPF-propagation involves a network propagation
#'        process based on random walk with restart (RWR) on a given gene-gene
#'        interaction network to obtain a distribution for each node (gene),
#'        which captures its relevance to all other genes in the network.
#'
#' @param x    A expression count matrix. The rows correspond to genes and
#'            the columns correspond to cells.
#' @param network  A adjacency matrix contation gene-gene interaction network.
#' @param gammma  A number between 0 and 1 (default: 0.5). \code{gamma} is the trade-off between prior
#'                information and network diffusion, governing the distance that a signal
#'                is allowed to diffuse through the network during smoothing.
#'                The specific value of \code{gamma} has little effect on the results of network
#'                propagation over a sizable range.
#' @param nThreads The number of cores to use. Default is 1.
#' @return A network-smoothed gene expression matrix.
#'
#' @import foreach
#'
#' @export
smoothRWR<- function(x,network="context",gamma=0.5,nThreads=1){
  if (foreach::getDoParWorkers() > 1) {
    if (nThreads > 1 & foreach::getDoParWorkers() != nThreads) {
      message(paste("Parallel backend already registered and is inconsistent",
                    "with ncores."))
    }
    nThreads <- getDoParWorkers()
  } else {
    if (nThreads > 1) {
      cl.create <- TRUE
      cl <- parallel::makeCluster(nThreads, outfile = "")
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      })
    }
  }#end else

  ep.data <- matrix(rep(0,nrow(network)*dim(x)[2]),
                    nrow=nrow(network))
  rownames(ep.data) <- rownames(network)
  colnames(ep.data) <- colnames(x)

  gene.overlap <- intersect(rownames(ep.data),
                            rownames(x))
  ep.data[gene.overlap,] <- x[gene.overlap,]
  Anorm <- network/Matrix::rowSums(network)
  eye <- diag(dim(network)[1])
  AA <- Matrix::t(eye - gamma*Anorm)
  BB <- (1-gamma) * ep.data
  smooth.data <- solve(AA, BB)
  smooth.x <- x
  smooth.x[gene.overlap,]<-as.matrix(smooth.data[gene.overlap,])
  return(smooth.x)
}




#' @title  Imputing dropout values of scRNA-seq data.
#'
#' @description scNPF-propagation for imputing dropouts and correcting expression
#'         expression measurements.scNPF-propagation involves a network propagation
#'        process based on random walk with restart (RWR) on a given gene-gene
#'        interaction network to obtain a distribution for each node (gene),
#'        which captures its relevance to all other genes in the network.
#'
#' @param x    A expression count matrix. The rows correspond to genes and
#'            the columns correspond to cells.
#' @param network  A adjacency matrix contation gene-gene interaction network.
#'                  User can use priori mode or context mode. For priori mode,
#'                  users can use publicly available molecular networks. In this
#'                  package, we provided three human gene-gene interaction networks,
#'                  including String, HumanNet and an integrated network. For context
#'                  mode (default), a context-specific gene-gene network is constructed from
#'                  the scRNA-seq data by WGCNA package.
#' @param gammma  A number between 0 and 1 (default: 0.5). \code{gamma} is the trade-off between prior
#'                information and network diffusion, governing the distance that a signal
#'                is allowed to diffuse through the network during smoothing.
#'                The specific value of \code{gamma} has little effect on the results of network
#'                propagation over a sizable range.
#' @param qt.gene A numeric value between 0 and 1 (default: 0.4) indicating the
#'                top percent of expressed genes to be reserved for buliding a context-specific
#'                gene-gene network. Used only if \code{network} = "context".
#' @param qt.cell A numeric value between 0 and 1 (default: 0.5) indicating the
#'                top percent of expressed cells to be reserved for buliding a context-specific
#'                gene-gene network. Used only if \code{network} = "context".
#' @param nThreads The number of cores to use. Default is 1.
#' @return A network-smoothed gene expression matrix.
#'
#' @examples
#' #Loading example scRNA-seq data.
#' load(system.file("data","yan.Rdata",package = "scNPF"))
#' #Testing with all genes
#' exp.data <- yan$data
#' #Or testing with paritial genes
#' #exp.data <- yan$data[1:2000,]
#'
#' ##For context mode
#' context.data <- scNPF.pro(x=exp.data, network="context",nThreads=8)
#' dim(context.data)
#' dim(exp.data)
#' context.data[1:5,1:3]
#' exp.data[1:5,1:3]
#'
#' ##For priori mode
#' ##Using String network
#' load(system.file("data","string.Rdata",package = "scNPF"))
#' string.data <- scNPF.pro(x=exp.data, network=string,nThreads=8)
#'
#' ##Using HumanNet network
#' load(system.file("data","humannet.Rdata",package = "scNPF"))
#' hm.data <- scNPF.pro(x=exp.data,network=humannet,nThreads=8)
#'
#' ##Using integrated network
#' load(system.file("data","integrated.Rdata",package = "scNPF"))
#' inter.data <- scNPF.pro(x=exp.data,network=INet,nThreads=8)
#'
#' @export
scNPF.pro <- function(x, network, gamma=0.5,qt.gene=0.4,qt.cell=0.5,nThreads=1) {
          if(!is.matrix(x)){
              stop( "'x' must be a expression count matrix")
            }
          if(!((is.numeric(gamma) && (gamma > 0 && gamma < 1)))){
              stop("'gamma' must be between 0 and 1")
           }
           if(!((is.numeric(qt.gene) && (qt.gene >= 0 && qt.gene <= 1)))){
             stop("'qt.gene' must be between 0 and 1")
           }
          if(!((is.numeric(qt.cell) && (qt.cell >= 0 && qt.cell <= 1)))){
             stop("'qt.cell' must be between 0 and 1")
          }
         if(is.character(network)){
            if(network == "context"){
                network <- context.mode(data=x,qt.gene=qt.gene,qt.cell=qt.cell,nThreads=nThreads)
            }else{
              stop("without this mode")
            }
        }
          if(!(is(network, 'matrix') || is(network, 'sparseMatrix'))){
              stop("'network' must be  a adjacency matrix" )
          }
          rownames(x) <- toupper(rownames(x))
          rownames(network) <- toupper(rownames(network) )
          colnames(network) <- toupper(colnames(network) )
          index <- match(rownames(network),rownames(x))
          network <- network[!is.na(index),!is.na(index)]
          nullrows <- Matrix::rowSums(network)==0
          network <- network[!nullrows,!nullrows]
          if(sum(Matrix::rowSums(network)==0)>0) {
            stop("gene interaction network cannot have zero rows/columns")
          }
          if(sum(Matrix::colSums(network)==0)>0) {
            stop("gene interaction network cannot have zero rows/columns")
          }
          x.smoothed <- smoothRWR(x=x, network=network, gamma=gamma,nThreads=nThreads)
          return(x.smoothed)
}



