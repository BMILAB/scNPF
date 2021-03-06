% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scNPF_fus.R
\name{scNPF.fus}
\alias{scNPF.fus}
\title{Similarity matrix Construction}
\usage{
scNPF.fus(data, K = 20, alpha = 0.5, T = 10)
}
\arguments{
\item{data}{A list stored network-smoothed gene expression
matrix with two length.}

\item{K}{Number of neighbors in K-nearest neighbors part of the algorithm.}

\item{T}{Number of iterations for the diffusion process.}

\item{alpah}{Variance for local model.}
}
\value{
A consensus similarity matrix.
}
\description{
constructs a sample-similarity network
     for each propagated expression matrix and then
     integrates these networks into a single cell-cell
     similarity network based on a nonlinear combination
     method.
}
\examples{
##Loading example data
load(system.file("data","yan.Rdata",package = "scNPF"))
exp.data <- yan$data[1:2000,]

##Smoothed gene expression based on priori mode
##Using String network.
load(system.file("data","string.Rdata",package = "scNPF"))
string.data <- scNPF.pro(x=exp.data, network=string)

##Smmothed gene expression based on context mode
context.data<- scNPF.pro(x=exp.data, network="context",qt.gene=0.9,qt.cell=0.9)

##Construction a cell-by-cell similarity matrix.
similarity <- scNPF.fus(data=list(string=string.data,context=context.data))

}
