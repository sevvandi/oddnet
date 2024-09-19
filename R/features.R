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
#' @return A network features object containing 20 graph-theoretic features.
#'
#' @examples
#' set.seed(1)
#' gr <- igraph::sample_gnp(100, 0.05)
#' compute_features(gr)
#'
#'
#' @importFrom stats median quantile sd
#' @export compute_features
compute_features <- function(gr,
                             attributes = FALSE,
                             attr_name = NULL,
                             fast = FALSE){
  # gr is an igraph object
  # triangles, degree and edges
  # compute_features_5
  triangles <- igraph::count_triangles(gr)
  triangles_50 <- median(triangles)
  triangles_99 <- quantile(triangles, probs = 0.99)
  degobj <- igraph::degree(gr)
  degree_50 <- median(degobj)
  degree_99 <- quantile(degobj, probs = 0.99)
  edges <- igraph::gsize(gr)
  density <- igraph::edge_density(gr)
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
    uniq_vals <- unique(igraph::vertex_attr(gr, vert_attr))
    num_set <- 1:length(uniq_vals)
    assortativity <-  igraph::assortativity.nominal(gr, types = num_set[as.factor(igraph::vertex_attr(gr, vert_attr))])
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
  clust_size_99 <- quantile(clust_obj, probs = 0.99)
  # centrality measures
  close_obj <- igraph::closeness(gr)
  closeness_50 <- median(close_obj, na.rm = TRUE)
  closeness_99 <- quantile(close_obj, na.rm = TRUE, probs = 0.99)
  close_overpoint8 <- sum(close_obj >= 0.8, na.rm = TRUE)/length(close_obj)
  between_obj <- igraph::betweenness(gr)
  between_50 <- median(between_obj, na.rm = TRUE)
  between_99 <- quantile(between_obj, na.rm = TRUE, probs = 0.99)
  pagerank_obj <- igraph::page_rank(gr)$vector
  pagerank_50 <- median(pagerank_obj)
  pagerank_99 <- quantile(pagerank_obj, probs = 0.99)
  cores_obj <- igraph::coreness(gr)
  cores_50 <- median(cores_obj)
  cores_99 <- quantile(cores_obj, probs = 0.99)
  hubs_obj <- igraph::hub_score(gr, scale = FALSE)
  hubs <- hubs_obj$value
  authority_obj <- igraph::authority_score(gr, scale = FALSE)
  authority <- authority_obj$value
  structure(list(
    num_nodes = num_nodes,
    triangles_99 = triangles_99,
    degree_99 = degree_99,
    edges = edges,
    density = density,
    clustering_coef = clust_coeff,
    assortativity = assortativity,
    mean_distance = mean_dist,
    diameter = diam,
    isolates = isolates,
    connectivity = conectivity,
    efficiency = efficiency,
    num_clusters = num_clust,
    cluster_size_99 = clust_size_99,
    close_overpoint8 = close_overpoint8,
    betweenness_99 = between_99,
    pagerank_99 = pagerank_99,
    cores_99 = cores_99,
    hubs_value = hubs,
    authority_value = authority,
    call = match.call()),
    class='networkfeatures')
}
