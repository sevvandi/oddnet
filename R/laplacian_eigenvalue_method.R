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

#' Laplacian Eigen Value method by Shenyang Huang, Yasmeen Hitti, Guillaume Rabusseau
#' and Reihaneh Rabbany from their KDD'20 paper Laplacian Change Point Detection for Dynamic Graphs
#'
#' @param matlist The matrix list, where each matrix is an adjacency matrix of the graph.
#' @param k The number of eigen values to connsider
#' @param short_win The length of the shorter windows
#' @param long_win The length of the longer windows
#' @param alpha The threshold to declare anomalies
#' @param from_file This is an additional parameter only if a file needs to be read
#' @importFrom utils read.csv
#'
#' @returns An object of class lad. LAD is a window based method. It considers short and a long
#' windows.  The lad object has anomalous scores when taking into account short and long windows
#' along with the identified anomalies for both short and long windows.
#'
#'
#' @examples
#' # We generate a series of networks and add an anomaly at 50th network.
#' set.seed(1)
#' networks <- list()
#' p.seq <- rep(0.05, 50)
#' p.seq[20] <- 0.2  # anomalous network at 20
#' for(i in 1:50){
#'   gr <- igraph::sample_gnp(100, p = p.seq[i])
#'   networks[[i]] <- igraph::as_adjacency_matrix(gr)
#' }
#' ladobj <- lad(networks, k = 6, short_win = 2, long_win = 4)
#' ladobj
#'
#'@references Huang, S., Hitti, Y., Rabusseau, G., & Rabbany, R. (2020). Laplacian Change
#'Point Detection for Dynamic Graphs. Proceedings of the ACM SIGKDD International Conference
#'on Knowledge Discovery and Data Mining, 349–358. https://doi.org/10.1145/3394486.3403077

#'
#' @export
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


