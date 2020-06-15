#' @title Plot Heatmap
#'
#' @description This function plots a heatmap figure form a distance or
#'   similarity matrix. The colors denote the sample types based on \code{label}.
#'
#' @param dists A distance or similarity matrix.
#' @param label True labels of sample.
#' @return This function returns a heatmap object invisibly.
#' @examples
#' library(gplots)
#' library(RColorBrewer)
#' load(system.file("data","yan.Rdata",package = "scNPF"))
#' data.dist <- as.matrix(dist(t(yan$data)))
#' #dev.off()
#' plotHeatmap(data.dist,yan$label)
#' @export
plotHeatmap <- function(dists,label){
  dists <- Turn_distance(dists)
  label <- as.character(label)
  index <- order(label)
  type <- label[index]
  y <- type
  for (i in 1:length(unique(type))){
    y[which(y==unique(type)[i])] <-i

  }
  y <- as.numeric(y)
  cols <- palette(brewer.pal(length(unique(y)), "Set1"))[y]
  heatmap.2(as.matrix(dists)[index,index],trace="none",
              Rowv=FALSE, Colv=FALSE,
              ColSideColors=cols, labRow=FALSE, labCol=FALSE,
              RowSideColors=cols,
              key.title=TRUE,dendrogram="none",key=TRUE,
              density.info="none",main=names(dists)[i],keysize = 1)

  dev.off()
}


