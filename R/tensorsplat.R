tensorsplat <- function(matlist, k = 5, alpha = 0.05 ){
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

  # for(i in 1:length(matlist)){
  #   numrows <- dim(matlist[[i]])[1]
  #   network_tensor[i, 1:numrows, 1:numrows] <- matlist[[i]] # as.matrix(matlist[[i]])
  # }

  cp_decomp <- rTensor::cp(network_tensor, num_components=3)
  temp_factors <- cp_decomp$U[[1]]
  lofout <- DDoutlier::LOF(temp_factors, k)
  alphath <- ceiling(length(lofout)*alpha)
  ord <- order(lofout, decreasing = TRUE)[1:alphath]
  ord
}


# tensorsplat2 <- function(matlist, k = 5 ){
#   nn <- dim(matlist[[1]])[1]
#   for(i in 1:length(matlist)){
#     numrows <- dim(matlist[[i]])[1]
#     nn <- max(nn, numrows)
#   }
#
#   indices <- c(length(matlist), nn, nn)
#   arr <- array(0, dim = indices)
#   for(i in 1:length(matlist)){
#     numrows <- dim(matlist[[i]])[1]
#     arr[i, , ] <- as.matrix(matlist[[i]])
#   }
#
#
#   cp_decomp <-multiway::parafac(arr, nfac = 3, nstart = 1)
#   temp_factors <- cp_decomp$U[[1]]
#   lofout <- DDoutlier::LOF(temp_factors, k)
#   lofout
# }
