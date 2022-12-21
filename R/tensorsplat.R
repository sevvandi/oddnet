#' Implements Danai Koutra's TensorSplat algorithm
#'
#' @param matlist The list of matrices where each matrix denotes the adjacency matrix of the network.
#' @param k The number of components in PARFAC tensor decomposition.
#' @param alpha The threshold to declare anomalies
#'
#' @returns Anomalous observations
#'
#' @examples
#' # We generate a series of networks and add an anomaly at 50th network.
#' set.seed(1)
#' networks <- list()
#' p.or.m.seq <- rep(0.05, 50)
#' p.or.m.seq[20] <- 0.2  # anomalous network at 20
#' for(i in 1:100){
#'   gr <- igraph::erdos.renyi.game(100, p.or.m = p.or.m.seq[i])
#'   networks[[i]] <- igraph::as_adjacency_matrix(gr)
#' }
#' tensobj <- tensorsplat(networks, k = 2)
#' tensobj  # anomalous networks
#'
#'
#' @references Koutra, D., Papalexakis, E. E., & Faloutsos, C. (2012).
#' TensorSplat: Spotting latent anomalies in time. Proceedings of the 2012 16th
#' Panhellenic Conference on Informatics, PCI 2012, 144–149.
#' https://doi.org/10.1109/PCi.2012.60
#'
#' @export
tensorsplat <- function(matlist, k = 2, alpha = 0.05 ){
  nn <- dim(matlist[[1]])[1]
  for(i in 1:length(matlist)){
    numrows <- dim(matlist[[i]])[1]
    nn <- max(nn, numrows)
  }

  indices <- c(length(matlist), nn, nn)
  arr <- array(0, dim = indices)
  for(i in 1:length(matlist)){
    numrows <- dim(matlist[[i]])[1]
    arr[i, 1:numrows, 1:numrows] <- as.matrix(matlist[[i]])
  }
  network_tensor <- rTensor::as.tensor(arr);

  cp_decomp <- rTensor::cp(network_tensor, num_components=3)
  temp_factors <- cp_decomp$U[[1]]
  lofout <- DDoutlier::LOF(temp_factors, k)
  alphath <- ceiling(length(lofout)*alpha)
  ord <- order(lofout, decreasing = TRUE)[1:alphath]
  ord
}



