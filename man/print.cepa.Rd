\name{print.cepa}
\alias{print.cepa}
\title{print the cepa object}
\description{
print the general description of a cepa object
}
\usage{
\method{print}{cepa}(x, ...)                     
}
\arguments{
  \item{x}{a \link{cepa} object}
  \item{...}{other arguments}
}
\details{
The function print the number of differential nodes in the pathway,
the number of all nodes in the pathway, the number of differential
genes in the pathway and the number of all genes in the pathway.

Also, the centrality measurement, the p-value are printed.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}
}
\examples{
data(PID.db)
data(gene.list)

i = 10
# get the 10th pathway in NCI catalogue
path = PID.db$NCI$pathList[[i]]
# get the interactions in this pathway
inter = PID.db$NCI$interactionList[PID.db$NCI$interactionList[, 1] \%in\% path, 2:3]
# ensure the interaction list is a matrix
inter = as.matrix(inter)
# find whether this pathway is significant
res = cepa(gene.list$dif, gene.list$bk, generate.pathway(inter),
           PID.db$NCI$mapping, "in.degree", iter = 500)
res
}
