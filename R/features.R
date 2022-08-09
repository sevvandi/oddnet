#' Computes features for each network.
#'
#' This function computes features for each network using graph theoretic constructs.
#'
#' @param gr The network or graph as an \code{igraph} object.
#' @param attributes If the network nodes/vertices have attributes, then \code{attributes = TRUE}.
#' @param attr_name The name of the node/vertex attribute. Only a single attribute can
#' be specified.
#' @param fast If set to \code{TRUE} will avoid computing time consuming features.
#'
#' @return A network features object.
#'
#' @examples
#' set.seed(1)
#' gr <- igraph::erdos.renyi.game(100, 0.05)
#' compute_features_4(gr)
#'
#'
#' @importFrom stats median quantile
#' @export compute_features_4
compute_features_4 <- function(gr, attributes = FALSE, attr_name = NULL, fast = FALSE){
  # gr is an igraph object
  # triangles, degree and edges

  triangles <- igraph::count_triangles(gr)
  triangles_50 <- median(triangles)
  triangles_95 <- quantile(triangles, probs = 0.95)
  degobj <- igraph::degree(gr)
  degree_50 <- median(degobj)
  degree_95 <- quantile(degobj, probs = 0.95)
  edges <- igraph::gsize(gr)
  num_nodes <- length(degobj)
  # clustering coefficient
  clust_coeff <- igraph::transitivity(gr)

  # assortativity
  if(attributes){
    vert_attr <- names(igraph::vertex_attr(gr))
    if(length(vert_attr) > 1){
      if(is.null(attr_name)){
        stop("Too many vertex attributes. Specify a single attribute using parameter attr_name for feature computations.")
      }else{
        vert_attr <- attr_name
      }
    }
    assortativity <-  igraph::assortativity.nominal(gr, types = c(1,2)[as.factor(igraph::vertex_attr(gr, vert_attr))])
  }else{
    assortativity <- igraph::assortativity_degree(gr)
  }
  # distance and isolates and connectivity
  mean_dist <- igraph::mean_distance(gr)
  diam <- igraph::diameter(gr)
  isolates <- sum(degobj == 0)/length(degobj)  # isolated nodes percentage
  if(!fast){
    conectivity <- igraph::cohesion(gr)
  }else{
    conectivity <- NA
  }

  # efficiency
  efficiency <- igraph::global_efficiency(gr)
  # number of clusters and cluster size
  clust <- igraph::clusters(gr)
  num_clust <- clust$no
  clust_obj <- clust$csize
  clust_size_50 <- median(clust_obj)
  clust_size_95 <- quantile(clust_obj, probs = 0.95)
  # centrality measures
  close_obj <- igraph::closeness(gr)
  closeness_50 <- median(close_obj, na.rm = TRUE)
  closeness_95 <- quantile(close_obj, na.rm = TRUE, probs = 0.95)
  between_obj <- igraph::betweenness(gr)
  between_50 <- median(between_obj, na.rm = TRUE)
  between_95 <- quantile(between_obj, na.rm = TRUE, probs = 0.95)
  pagerank_obj <- igraph::page_rank(gr)$vector
  pagerank_50 <- median(pagerank_obj)
  pagerank_95 <- quantile(pagerank_obj, probs = 0.95)
  cores_obj <- igraph::coreness(gr)
  cores_50 <- median(cores_obj)
  cores_95 <- quantile(cores_obj, probs = 0.95)
  structure(list(
    num_nodes = num_nodes,
    triangles_95 = triangles_95,
    degree_95 = degree_95,
    edges = edges,
    clustering_coef = clust_coeff,
    assortativity = assortativity,
    mean_distance = mean_dist,
    diameter = diam,
    isolates = isolates,
    connectivity = conectivity,
    efficiency = efficiency,
    num_clusters = num_clust,
    cluster_size_95 = clust_size_95,
    closeness_95 = closeness_95,
    betweenness_95 = between_95,
    pagerank_95 = pagerank_95,
    cores_95 = cores_95,
    call = match.call()),
    class='networkfeatures')
}

