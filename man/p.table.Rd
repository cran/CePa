\name{p.table}
\alias{p.table}
\title{Table of p-values of pathways}
\description{
Table of p-values of pathways
}
\usage{
p.table(x)
}
\arguments{
  \item{x}{\link{cepa.all} object}
}
\value{
A data matrix where rows represent pathways and columns represent centralities.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}
}
\examples{
data(PID.db)
data(gene.list)
res = cepa.all(gene.list$dif, gene.list$bk, PID.db$NCI$pathList[1:2],
               PID.db$NCI$interactionList, PID.db$NCI$mapping, iter = 500)
p.table(res)
}
