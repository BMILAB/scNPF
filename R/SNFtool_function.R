#' @title Affinity matrix calculation
#'
#' @description Computes affinity matrix from a generic distance matrix
#'
#'
#' @param Diff Distance matrix
#' @param K Number of nearest neighbors
#' @param sigma Variance for local model
#' @return Returns an affinity matrix that represents the neighborhood graph of
#'        the data points.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @references B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B
#'         Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and
#'         effective method to aggregate multiple data types on a genome wide scale.
#'         Nature Methods. Online. Jan 26, 2014
#' @examples
#'
#'
#' ## First, set all the parameters:
#' K = 20; ##number of neighbors, must be greater than 1. usually (10~30)
#' alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
#' T = 20; ###Number of Iterations, usually (10~50)
#'
#' ## Data1 is of size n x d_1,
#' ## where n is the number of patients, d_1 is the number of genes,
#' ## Data2 is of size n x d_2,
#' ## where n is the number of patients, d_2 is the number of methylation
#' data(Data1)
#' data(Data2)
#'
#' ## Calculate distance matrices(here we calculate Euclidean Distance,
#' ## you can use other distance, e.g. correlation)
#' Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
#' Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))
#'
#' ## Next, construct similarity graphs
#' W1 = affinityMatrix(Dist1, K, alpha)
#' W2 = affinityMatrix(Dist2, K, alpha)
#'
#' @export
affinityMatrix <- function(Diff,K=20,sigma=0.5) {
  ###This function constructs similarity networks.
  N <- nrow(Diff)
  Diff <- (Diff + t(Diff)) / 2
  diag(Diff) <- 0;
  sortedColumns <- as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }
  means <- apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
  avg <- function(x,y) ((x+y)/2)
  Sig <- outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
  densities <- dnorm(Diff,0,sigma*Sig,log = FALSE)
  W <- (densities + t(densities)) / 2
  return(W)
}

#' @title Similarity Network Fusion
#'
#' @description Similarity Network Fusion takes multiple views of a network and fuses them
#' together to construct an overall status matrix. The input to our algorithm
#' can be feature vectors, pairwise distances, or pairwise similarities. The
#' learned status matrix can then be used for retrieval, clustering, and
#' classification.
#'
#'
#' @param Wall List of matrices. Each element of the list is a square,
#' symmetric matrix that shows affinities of the data points from a certain
#' view.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @param layer_bias Numerical vector of length length(Wall) containing the
#' relative weights for each layer.
#' @param parallel Should the algorithm attempt to run in parallel? A parallel
#' backend needs to be set up first.
#' @param auto_stop Should the algorithm stop early if it believes it has
#' converged?
#' @return W is the overall status matrix derived
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @references B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B
#' Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and
#' effective method to aggregate multiple data types on a genome wide scale.
#' Nature Methods. Online. Jan 26, 2014
#'
#' Concise description can be found here:
#' http://compbio.cs.toronto.edu/SNF/SNF/Software.html
#' @examples
#'
#'
#' ## First, set all the parameters:
#' K = 20;		# number of neighbors, usually (10~30)
#' alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
#' T = 20; 	# Number of Iterations, usually (10~20)
#'
#' ## Data1 is of size n x d_1,
#' ## where n is the number of patients, d_1 is the number of genes,
#' ## Data2 is of size n x d_2,
#' ## where n is the number of patients, d_2 is the number of methylation
#' data(Data1)
#' data(Data2)
#'
#' ## Here, the simulation data (SNFdata) has two data types. They are complementary to each other.
#' ## And two data types have the same number of points.
#' ## The first half data belongs to the first cluster; the rest belongs to the second cluster.
#' truelabel = c(matrix(1,100,1),matrix(2,100,1)); ## the ground truth of the simulated data
#'
#' ## Calculate distance matrices
#' ## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
#'
#' ## If the data are all continuous values, we recommend the users to perform
#' ## standard normalization before using SNF,
#' ## though it is optional depending on the data the users want to use.
#' # Data1 = standardNormalization(Data1);
#' # Data2 = standardNormalization(Data2);
#'
#'
#'
#' ## Calculate the pair-wise distance;
#' ## If the data is continuous, we recommend to use the function "dist2" as follows
#' Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
#' Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
#'
#' ## next, construct similarity graphs
#' W1 = affinityMatrix(Dist1, K, alpha)
#' W2 = affinityMatrix(Dist2, K, alpha)
#'
#' ## These similarity graphs have complementary information about clusters.
#' displayClusters(W1,truelabel);
#' displayClusters(W2,truelabel);
#'
#' ## next, we fuse all the graphs
#' ## then the overall matrix can be computed by similarity network fusion(SNF):
#' W = SNF(list(W1,W2), K, T)
#'
#' ## With this unified graph W of size n x n,
#' ## you can do either spectral clustering or Kernel NMF.
#' ## If you need help with further clustering, please let us know.
#'
#' ## You can display clusters in the data by the following function
#' ## where C is the number of clusters.
#' C = 2 								# number of clusters
#' group = spectralClustering(W,C); 	# the final subtypes information
#' displayClusters(W, group)
#'
#' ## You can get cluster labels for each data point by spectral clustering
#' labels = spectralClustering(W, C)
#'
#' plot(Data1, col=labels, main='Data type 1')
#' plot(Data2, col=labels, main='Data type 2')
#'
#' @export
SNF <- function(Wall,K=20,t=20,layer_bias = rep.int(1, length(Wall)),parallel=FALSE, auto_stop=FALSE) {

  ###This function is the main function of our software. The inputs are as follows:
  # Wall : List of affinity matrices
  # K : number of neighbors
  # t : number of iterations for fusion
  # layer_bias : Numerical vector of length length(Wall) containing the relative weights for each layer.
  # parallel : Should the algorithm attempt to run in parallel? A parallel backend needs to be set up first.
  # auto_stop : Should the algorithm stop early if it believes it has converged?

  ###The output is a unified similarity graph. It contains both complementary information and common structures from all individual network.
  ###You can do various applications on this graph, such as clustering(subtyping), classification, prediction.

  LW <- length(Wall)
  normalize <- function(X) X / rowSums(X)
  # makes elements other than largest K zero

  newW <- vector("list", LW)
  nextW <- vector("list", LW)
  ###First, normalize different networks to avoid scale problems.
  for( i in 1: LW){
    Wall[[i]] <- normalize(Wall[[i]]);
    Wall[[i]] <- (Wall[[i]]+t(Wall[[i]]))/2;
  }

  ### Calculate the local transition matrix.
  for( i in 1: LW){
    newW[[i]] <- (.dominateset(Wall[[i]],K))
  }

  # Set up convergence detection
  devs <- NULL
  converged <- FALSE

  # perform the diffusion for t iterations
  for (i in 1:t) {
    nextW <- plyr::llply(1:LW, function(j){
      sumWJ <- matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      adjusted_bias <- layer_bias/sum(layer_bias[-j])
      for (k in setdiff(1:LW, j)) {
        sumWJ <- sumWJ + (Wall[[k]]*adjusted_bias[k])
      }
      return(newW[[j]] %*% (sumWJ/(LW - 1)) %*% t(newW[[j]]))
    }, .parallel=parallel)
    ###Normalize each new obtained networks.
    for(j in 1 : LW){

      Wall[[j]] <- nextW[[j]] + diag(nrow(Wall[[j]]));
      Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2;
    }


    # Check for convergence
    wallarray <- array(unlist(Wall),dim=c(dim(Wall[[1]]),length(Wall)))
    means <- base::rowMeans(wallarray, dims=2)
    dev <- mean(abs(wallarray-array(means,dim=dim(wallarray))) / mean(means))
    devs <- c(devs,dev)
    d_devs <- (devs-lag(devs))
    d2_devs <- (d_devs-lag(d_devs))
    if(!is.na(d2_devs[length(d2_devs)]) &
       abs(d_devs[length(d_devs)])<0.01 &
       abs(d2_devs[length(d2_devs)])<0.01){
      converged <- TRUE
    }
    if(auto_stop & converged){
      break
    }
  }

  # construct the combined affinity matrix by summing diffused matrices
  W <- matrix(0,nrow(Wall[[1]]), ncol(Wall[[1]]))
  for(i in 1:LW){
    W <- W + Wall[[i]]
  }
  W <- W/LW;
  W <- normalize(W);
  # ensure affinity matrix is symmetrical
  W <- (W + t(W)+diag(nrow(W))) / 2;

  return(W)
}


#' @export
.dominateset <- function(xx,KK=20) {
  ###This function outputs the top KK neighbors.

  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);

  }
  return(normalize(A))
}

#' @export
Turn_distance <- function(W){
  normalize <- function(X) X/rowSums(X)
  diag(W) = median(as.vector(W))
  W = normalize(W)
  W = W + t(W)
  return(W)
}
