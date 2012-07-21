\name{cepa.all.parallel}
\alias{cepa.all.parallel}
\title{use CePa package by parallel computing}
\description{
use CePa package by parallel computing
}
\usage{
cepa.all.parallel(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    glevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2)  
}
\arguments{
  \item{dif}{differential gene list}
  \item{bk}{background gene list}
  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \link{sampleLabel} object identify the design of the microarray experiment}
  \item{pc}{a \link{pathway.catalogue} object storing information of pathways}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names}
  \item{glevel}{gene level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs"}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank"}
  \item{iter}{number of simulations}
  \item{ncores}{number of cores for parallel computing}
}
\details{
The function divide the pathway list into several parts and each part is sent to a core for 
parallel computing.

The package for parallel computing is \code{snow}.

Note: there may be warnings saying connections not closed. I don't know why this
happens. However it does not matter unless you have obsessive compulsive disorder.
}
\value{
    A \code{\link{cepa.all}} class object
}
\references{
Gu Z, Liu J, Cao K, Zhang J, Wang J. Centrality-based pathway enrichment: a systematic 
approach for finding significant pathways dominated by key genes. BMC Syst Biol. 2012 Jun 6;6(1):56.
}
\author{Zuguang Gu <jokergoo@gmail.com>}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{

data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all.parallel(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI, ncores = 4)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI, ncores = 4)
}
}
