#' @title Building context mode
#'
#' @description a context-specific gene-gene network is constructed from
#'                  the scRNA-seq data by WGCNA package.
#'
#' @param data    A expression count matrix. The rows correspond to genes and
#'            the columns correspond to cells.
#'
#' @param qt.gene A numeric value between 0 and 1 (default: 0.4) indicating the
#'                top percent of expressed genes to be reserved for buliding a context-specific
#'                gene-gene network. Used only if \code{network} = "context".
#' @param qt.cell A numeric value between 0 and 1 (default: 0.5) indicating the
#'                top percent of expressed cells to be reserved for buliding a context-specific
#'                gene-gene network. Used only if \code{network} = "context".
#' @param nThreads The number of cores to use. Default is 1.
#'
#' @examples
#' #Loading example scRNA-seq data.
#' load(system.file("data","yan.Rdata",package = "scNPF"))
#' exp.data <- yan$data
#' dim(exp.data)
#' label <- yan$label
#'
#' #Parallel
#' context.network <- context.mode(data=exp.data,nThreads=8)
#' dim(context.network)
#' class(context.network)
#' @return A context-specific gene-gene network.
#'
#' @import foreach
#'
#' @export
context.mode<- function(data,qt.gene=0.4,qt.cell=0.5,nThreads=1){
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
  data <- log2(data+2)
  rownames(data) <- toupper(rownames(data))
  qt.gene <- 1- qt.gene
  qt.cell <- 1-qt.cell
  UMI  <- rowSums(data)
  cutoff <- quantile(UMI,qt.gene)
  index <- which(UMI>cutoff)
  data <- data[index,]
  count <- colSums(data)
  cutoff <- quantile(count,qt.cell)
  index <- which(count>cutoff)
  data <- data[,index]
  type <- "unsigned"
  corType <- "pearson"
  m.mad <- apply(data,1,mad)
  data <- data[which(m.mad>0),]
  data <- as.data.frame(t(data))
  gsg <- goodSamplesGenes(data, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(data)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(data)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    data <- data[gsg$goodSamples, gsg$goodGenes]
  }

  nGenes <- ncol(data)
  nSamples <- nrow(data)
  powers <- c(c(1:10), seq(from = 12, to=30, by=2))
  sft <- WGCNA::pickSoftThreshold(data, powerVector=powers,
                                  networkType=type, verbose=5)
  power <- sft$powerEstimate
  if (is.na(power)){
    power <- ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                    ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                           ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                  ifelse(type == "unsigned", 6, 12))
                    )
    )
  }
  TOM <- WGCNA::TOMsimilarityFromExpr(data, power=power, corType=corType, networkType=type,nThreads = nThreads)
  TOM <- as.matrix(TOM)
  rownames(TOM) <- colnames(data)
  colnames(TOM) <- colnames(data)

  smm <- TOM
  genename <- rownames(smm)
  gg <- igraph::graph_from_adjacency_matrix(smm, mode="undirected",weighted=TRUE)
  adj_matrix_d  <- igraph::as_adjacency_matrix(gg,attr="weight")
  rownames(adj_matrix_d) <- genename
  colnames(adj_matrix_d) <- genename
  adj_matrix_d<- adj_matrix_d[!duplicated(genename), !duplicated(genename)]
  nullrows <- Matrix::rowSums(adj_matrix_d )==0
  adj_matrix_d   <- adj_matrix_d [!nullrows,!nullrows]
  network <- adj_matrix_d
  return(network)
}


