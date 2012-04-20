\name{get.cepa}
\alias{get.cepa}
\title{get single cepa object from cepa.all object}
\description{
get single cepa object from cepa.all object
}
\usage{
get.cepa(x, id = NULL, cen = 1)    
}
\arguments{
  \item{x}{\link{cepa.all} object}
  \item{id}{index or the name of the pathway}
  \item{cen}{index or the name of the centrality}
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}, \code{\link{cepa.all}}
}
\examples{
data(PID.db)
data(gene.list)
res = cepa.all(gene.list$dif, gene.list$bk, PID.db$NCI$pathList[1:2],
               PID.db$NCI$interactionList, PID.db$NCI$mapping, iter = 500)
get.cepa(res, 2, "in.degree")
}
