\name{eccentricity}
\alias{eccentricity}
\title{Calculate eccentricity centrality}
\description{
Calculate eccentricity centrality
}
\usage{
eccentricity(graph, mode=c("all", "in", "out"))
}
\arguments{
  \item{graph}{\link[igraph]{igraph} object}
  \item{mode}{mode of the centrality}
}
\details{
The eccentricity identifies whether the node is in the center of the network.

The eccentricity is defined as 1/max(d(v, w)). where d(w, v) is the
 length of the shortest path from node w to node v.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}, \code{\link{spread}}, \code{\link{reach}}, \code{\link{radiality}}
}
\examples{
require(igraph)
pathway = barabasi.game(200)
eccentricity(pathway)
}
