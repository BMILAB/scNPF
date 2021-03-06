% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scNPF_pro.R
\name{smoothRWR}
\alias{smoothRWR}
\title{Imputing dropout values of scRNA-seq data.}
\usage{
smoothRWR(x, network = "context", gamma = 0.5, nThreads = 1)
}
\arguments{
\item{x}{A expression count matrix. The rows correspond to genes and
the columns correspond to cells.}

\item{network}{A adjacency matrix contation gene-gene interaction network.}

\item{nThreads}{The number of cores to use. Default is 1.}

\item{gammma}{A number between 0 and 1 (default: 0.5). \code{gamma} is the trade-off between prior
information and network diffusion, governing the distance that a signal
is allowed to diffuse through the network during smoothing.
The specific value of \code{gamma} has little effect on the results of network
propagation over a sizable range.}
}
\value{
A network-smoothed gene expression matrix.
}
\description{
scNPF-propagation for imputing dropouts and correcting expression
        expression measurements.scNPF-propagation involves a network propagation
       process based on random walk with restart (RWR) on a given gene-gene
       interaction network to obtain a distribution for each node (gene),
       which captures its relevance to all other genes in the network.
}
