\name{spread}
\alias{spread}
\title{
Calculate radiality centrality
}
\description{
Calculate radiality centrality
}
\usage{
spread(graph, mode = c("all", "in", "out"),
    weights = E(graph)$weight, f = function(x) 1/x)
}
\arguments{

  \item{graph}{an \code{\link[igraph:aaa-igraph-package]{igraph}} object}
  \item{mode}{mode of the centrality}
  \item{weights}{If edges in the graph have weight, then by default, the weight is used to calculate the length of the shortest path. Set it to NULL to supress the weight}
  \item{f}{function for the weaken rate}

}
\details{
The spread centrality measures how wide the node can send or receive the information in the network.
Like the water wave, the effect would be weakened with the increase of the distance to other nodes.

If the weaken function is defined as \code{1/x}, then the spread centrality is calculated as
\code{sum(1/d(w, v))} where \code{d(w, v)} is the length of the shortest path of node \code{w} and node \code{v}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\seealso{
\code{\link{reach}}, \code{\link{radiality}}
}
\examples{
require(igraph)
pathway = barabasi.game(200)
spread(pathway)
}
