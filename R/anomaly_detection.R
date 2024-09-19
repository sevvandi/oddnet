#' Identifies anomalous networks from a series of temporal networks.
#'
#' This function identifies anomalous networks from a series of temporal networks. It uses graph
#' theoretic features to transform networks to a feature space. This function has parameters for
#' feature computation, scaling, robust PCA and anomaly detection procedures.   ADD MORE DESCRIPTION.
#'
#' @param networks The input series of temporal networks given in a list with each network denoted
#' by its adjacency matrix.
#' @param alpha An anomaly detection parameter. The level of significance for the anomaly detection
#' algorithm lookout. Default is 0.05.
#' @param dd A robust PCA parameter. The number of reduced dimensions in robust PCA. Default is 2.
#' @param trim A scaling parameter. The percentage used to compute trimmed mean and trimmed standard
#' deviation. Default is 0.5 percent.
#' @param na_action The action for NA valued features.
#' @param vert_attr A feature computation parameter. If \code{TRUE} the network nodes/vertices have
#' attributes.
#' @param attr_name A feature computation parameter. The name of the network vertex attribute. Only a
#' single attribute can be specified.
#' @param attr_mat A feature computation parameter. If network nodes/vertices have attributes, the list of attribute
#'  matrices for each network can be given using this feature.
#' @param fast If set to \code{TRUE} will avoid computing time consuming features.
#' @param feature_subset A feature computation parameter. A subset of features to be used for anomaly detection.
#'
#' @return Object imported from lookout.
#' @seealso [lookout::lookout()]
#'
#' @examples
#' # We generate a series of networks and add an anomaly at 50th network.
#' set.seed(1)
#' networks <- list()
#' p.seq <- rep(0.1, 50)
#' p.seq[20] <- 0.3  # anomalous network at 20
#' for(i in 1:50){
#'   gr <- igraph::sample_gnp(50, p = p.seq[i])
#'   networks[[i]] <- igraph::as_adjacency_matrix(gr)
#' }
#' anomalous_networks(networks, fast = TRUE)
#'
#'
#' @importFrom dplyr '%>%'
#' @importFrom rlang .data
#'
#' @export anomalous_networks
anomalous_networks <- function(networks,
                               alpha = 0.05,
                               dd = 2,
                               trim = 0.005,
                               na_action = NULL,
                               vert_attr = FALSE,
                               attr_name = NULL,
                               attr_mat = NULL,
                               fast = FALSE,
                               feature_subset = NULL){


  num_networks <- length(networks)
  num_feat <- 20 # 26 # 26 for compute_features_6
  features <- matrix(0, nrow = num_networks, ncol = num_feat)

  # Compute network features
  for(i in 1:num_networks){
    gr <- igraph::graph_from_adjacency_matrix(networks[[i]])
    if(vert_attr){  # if the vertices have attributes
      gr <- igraph::set_vertex_attr(gr, attr_name, value = attr_mat[[i]])
      tt <- compute_features(gr, attributes = TRUE, attr_name = attr_name)
    }else{
      tt <- compute_features(gr, attributes = FALSE, attr_name = NULL, fast)
    }
    features[i, ] <- unlist(tt[1:num_feat])
  }
  colnames(features) <- names(tt)[1:num_feat]

  # If feature subset is given take only those features
  if(!is.null(feature_subset)){
    features <- features[ ,feature_subset]
  }

  if(!is.null(na_action)){
    features[is.na(features)] <- 0
  }
  # Remove features with all NA or NaN values
  na_sum <- apply(features, 2, function(x) sum(is.na(x)))
  inds <- which(na_sum == num_networks)
  if(length(inds) > 0){
    features <- features[ ,-inds]
  }
  # Remove features with constant values
  inds <- which(apply(features, 2, sd)==0)
  if(length(inds) > 0){
    features <- features[ ,-inds]
  }

  # Make a long data frame
  dfmerge <- tibble::as_tibble(features) %>%
    dplyr::mutate(t = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -t, values_to = "value", names_to = "feature") %>%
    tsibble::as_tsibble(index = t, key = .data$feature)

  # Fit ARIMA models
  fit <- dfmerge %>%
    fabletools::model(arima = fable::ARIMA(value))

  # Save residuals
  dfresiduals <- fabletools::augment(fit) %>%
    dplyr::select(t, .data$feature, .data$.innov) %>%
    tidyr::pivot_wider(names_from = .data$feature, values_from = .data$.innov) %>%
    tsibble::as_tibble() %>%
    dplyr::select(-t)

  dfresiduals2 <- scale_trimmed(dfresiduals, trim = trim)

  # PCA on all residuals
  pca <- pcaPP::PCAproj(dfresiduals2, scale = NULL, center = NULL, k = dd)
  dfpca <- tsibble::as_tibble(pca$scores) %>%
    dplyr::mutate(t = dplyr::row_number()) %>%
    tsibble::as_tsibble(index = t)

  # Find outliers in PCs using lookout
  lookobj <- lookout::lookout(dfpca[ ,1:dd], alpha = alpha)
  lookobj
}
