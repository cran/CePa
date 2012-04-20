\name{print.cepa.all}
\alias{print.cepa.all}
\title{print the cepa.all object}
\description{
print the cepa.all object
}
\usage{
\method{print}{cepa.all}(x, ...)                        
}
\arguments{
  \item{x}{\link{cepa.all} object}
  \item{...}{other arguments}
}
\details{
The function print the number of significant pathways under various
centrality measures at p-value < 0.01.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
data(PID.db)
data(gene.list)
# run the first 10 pathways as an example
res = cepa.all(gene.list$dif, gene.list$bk, PID.db$NCI$pathList[1:2],
               PID.db$NCI$interactionList, PID.db$NCI$mapping, iter = 500)
res

}
