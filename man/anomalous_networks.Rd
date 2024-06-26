% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anomaly_detection.R
\name{anomalous_networks}
\alias{anomalous_networks}
\title{Identifies anomalous networks from a series of temporal networks.}
\usage{
anomalous_networks(
  networks,
  alpha = 0.05,
  dd = 2,
  trim = 0.005,
  na_action = NULL,
  vert_attr = FALSE,
  attr_name = NULL,
  attr_mat = NULL,
  fast = FALSE,
  feature_subset = NULL
)
}
\arguments{
\item{networks}{The input series of temporal networks given in a list with each network denoted
by its adjacency matrix.}

\item{alpha}{An anomaly detection parameter. The level of significance for the anomaly detection
algorithm lookout. Default is 0.05.}

\item{dd}{A robust PCA parameter. The number of reduced dimensions in robust PCA. Default is 2.}

\item{trim}{A scaling parameter. The percentage used to compute trimmed mean and trimmed standard
deviation. Default is 0.5 percent.}

\item{na_action}{The action for NA valued features.}

\item{vert_attr}{A feature computation parameter. If \code{TRUE} the network nodes/vertices have
attributes.}

\item{attr_name}{A feature computation parameter. The name of the network vertex attribute. Only a
single attribute can be specified.}

\item{attr_mat}{A feature computation parameter. If network nodes/vertices have attributes, the list of attribute
matrices for each network can be given using this feature.}

\item{fast}{If set to \code{TRUE} will avoid computing time consuming features.}

\item{feature_subset}{A feature computation parameter. A subset of features to be used for anomaly detection.}
}
\value{
Object imported from lookout.
}
\description{
This function identifies anomalous networks from a series of temporal networks. It uses graph
theoretic features to transform networks to a feature space. This function has parameters for
feature computation, scaling, robust PCA and anomaly detection procedures.   ADD MORE DESCRIPTION.
}
\examples{
# We generate a series of networks and add an anomaly at 50th network.
set.seed(1)
networks <- list()
p.or.m.seq <- rep(0.1, 50)
p.or.m.seq[20] <- 0.3  # anomalous network at 20
for(i in 1:50){
  gr <- igraph::erdos.renyi.game(50, p.or.m = p.or.m.seq[i])
  networks[[i]] <- igraph::as_adjacency_matrix(gr)
}
anomalous_networks(networks, fast = TRUE)


}
\seealso{
[lookout::lookout()]
}
