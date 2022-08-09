context_matrix <- function(matlist, k = NULL){
  len <- length(matlist)
  if(is.null(k)){
    k <- min(10, ceiling(dim(matlist[[1]])[1]/2))
  }
  context <- matrix(0, nrow = k, ncol = len)
  for(ii in 1:len){
    adjacency <- matlist[[ii]]
    gr <- igraph::graph_from_adjacency_matrix(adjacency)
    # diag_degree <- diag(degree(gr))
    # laplacian <- diag_degree - adjacency
    laplacian <- igraph::laplacian_matrix(gr)
    svdout <- svd(laplacian, nu = k, nv = k)
    singvals <- svdout$d
    context[ ,ii] <- singvals[1:k]/(sqrt(sum(singvals[1:k]^2)))
  }
  context
}


lad_scores <- function(context, win_width ){
  timepoints <- dim(context)[2]
  num_times <- timepoints - win_width
  st <- 1
  en <- win_width
  scores <- rep(0, num_times)
  for(ii in 1:num_times){
    window <- context[ ,st:en]
    svdobj <- svd(window, nu = 1, nv = 1)
    norm_vec <- svdobj$u[ ,1]
    curr_vec <- context[ ,(en+1)]
    scores[ii] <- 1 - abs(sum(norm_vec*curr_vec)/sqrt(sum(norm_vec^2)*sum(curr_vec^2)))
    st <- st + 1
    en <- en + 1
  }
  scores
}


#' @importFrom utils read.csv
lad <- function(matlist, k = NULL, short_win, long_win, alpha = 0.05, from_file = NULL){
  # compute the context matrix
  if(!is.null(from_file)){
    filedat <- read.csv(from_file, header = FALSE)
    context <- t(filedat)
  }else{
    context <- context_matrix(matlist, k)
  }

  # computer short term scores
  anomaly_scores_short <- lad_scores(context, short_win)
  anomaly_scores_short <- c(rep(0, short_win), anomaly_scores_short)
  anomaly_scores_short <- c(0, diff(anomaly_scores_short))
  anomaly_scores_short[anomaly_scores_short < 0 ] <- 0

  anomaly_scores_long <- lad_scores(context, long_win)
  anomaly_scores_long <- c(rep(0, long_win), anomaly_scores_long)
  anomaly_scores_long <- c(0, diff(anomaly_scores_long))
  anomaly_scores_long[anomaly_scores_long < 0] <- 0

  alphath_pos <- ceiling(length(matlist)/20)
  anomalies_short <- (order(anomaly_scores_short, decreasing = TRUE))[1:alphath_pos]
  anomalies_long <- (order(anomaly_scores_long, decreasing = TRUE))[1:alphath_pos]

  structure(list(
    short_scores = anomaly_scores_short,
    long_scores = anomaly_scores_long,
    short_anomalies = anomalies_short,
    long_anomalies = anomalies_long,
    call = match.call()
  ), class='lad')

}
