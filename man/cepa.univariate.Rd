\name{cepa.univariate}
\alias{cepa.univariate}
\title{Access the significance of a pathway}
\description{
Using network centrality as weight to access the significance of a pathway.
Pathway nodes are taken as network units instead of genes.
}
\usage{
cepa.univariate(mat, label, pc, pathway = NULL, id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
    iter = 1000, glevel = "tvalue_abs", plevel = "mean",
    node.level.from.expr = NULL, node.level.t.value = NULL,
    r.node.level.from.expr = NULL)        
}
\arguments{
  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \link{sampleLabel} object identify the design of the microarray experiment}
  \item{pc}{a pathway.catalogue object storing information of pathways}
  \item{pathway}{\link[igraph]{igraph} object or edge list}
  \item{id}{identify the number of the pathway in the catalogue}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names}
  \item{glevel}{gene level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs"}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank"}
  \item{node.level.from.expr}{for simplicity of computing}
  \item{node.level.t.value}{for simplicity of computing}
  \item{r.node.level.from.expr}{for simplicity of computing}
  \item{iter}{number of simulations}
}
\value{
A \code{cepa} class object

}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{

data(PID.db)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.univariate(mat = eset, label = label, pc = PID.db$NCI, id = 2)
}
}