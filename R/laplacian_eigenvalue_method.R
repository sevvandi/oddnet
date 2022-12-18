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
#' @examples
#' \dontrun{
#' library(tnet)
#' library(lubridate)
#' data(tnet)
#' dat <- tnet::OnlineSocialNetwork.n1899.lnet
#' dat$t <-  as_datetime(dat$t)
#' dat$date <- as_date(dat$t)
#' length(which(dat$i == dat$j))
#' rminds <- which(dat$i == dat$j)
#' dat2 <- dat[-rminds, ]
#' dat3 <- dat2[ ,c("i", "j", "date")]
#' len <-  length(unique(dat3$date))
#' unique_dates <- unique(dat3$date)
#' num_networks <- length(unique_dates)
#' vset <- unique(sort(c(dat3[ ,1], dat3[ ,2])))
#' matlist <- list()
#' for(i in 1:20){
#'   nn <- unique_dates[i]
#'   inds <- which( dat3$date == nn )
#'   datwin <- dat3[inds, 1:2]
#'   gr <- graph_from_data_frame(datwin, vertices = vset)
#'   admat <- as_adjacency_matrix(gr)
#'   matlist[[i]] <- admat
#'}
#' ladobj <- lad(matlist, k = 6,  short_win = 7, long_win = 14)
#'}
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


