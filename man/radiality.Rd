\name{radiality}
\alias{radiality}
\title{Calculate radiality centrality}
\description{
Calculate radiality centrality
}
\usage{
radiality(graph, mode=c("all", "in", "out"))
}
\arguments{
  \item{graph}{\link[igraph]{igraph} object}
  \item{mode}{mode of the centrality}
}
\details{
The radiality is defined as sum(dG + 1 - d(v, w))/(n - 1). where d(w, v) is the
 length of the shortest path from node w to node v, dG is the diameter of the network,
 n is the size of the network.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}, \code{\link{spread}}, \code{\link{reach}}, \code{\link{radiality}}
}
\examples{
require(igraph)
pathway = barabasi.game(200)
radiality(pathway)
}
