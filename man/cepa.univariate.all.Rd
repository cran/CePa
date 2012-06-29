\name{cepa.univariate.all}
\alias{cepa.univariate.all}
\title{Apply centrality-extented GSA on a list of pathways}
\description{
Apply centrality-extented GSA on a list of pathways
}
\usage{
cepa.univariate.all(mat, label, pc, cen = default.centralities,
           cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
           glevel = "tvalue_abs", plevel = "mean", iter = 1000)    
}
\arguments{
  \item{mat}{expression matrix in which rows are genes and columns are samples}
  \item{label}{a \link{sampleLabel} object identify the design of the microarray experiment}
  \item{pc}{a pathway.catalogue object storing information of pathways}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names}
  \item{glevel}{gene level transformation, should be one of "tvalue", "tvalue_sq", "tvalue_abs"}
  \item{plevel}{pathway level transformation, should be one of "max", "min", "median", "sum", "mean", "rank"}
  \item{iter}{number of simulations}
}
\details{
The traditional gene-set analysis (GSA) to find significant pathways 
uses the whole expression matrix. GSA methods are implemented via either a univariate 
or a multivariate procedure. In univariate analysis, gene level statistics are 
initially calculated from fold changes or statistical tests (e.g., t-test). 
These statistics are then combined into a pathway level statistic by summation or 
averaging. Multivariate analysis considers the correlations between genes in the 
pathway and calculates the pathway level statistic directly from the expression 
value matrix using Hotelling's T^2 test or MANOVA models. The function implement
univariate procedure of GSA with network centralities.

If the centrality measurement is set as a string, only pre-defined "equal.weight",
"in.degree", "out.degree", "degree", "betweenness", "in.reach", "out.reach",
"reach", "in.spread", "out.spread" and "spread" are allowed. More centrality
measurements can be used by setting it as a function (such as closeness,
cluster coefficient). In the cepa.all function, we recommand users choose
at least two centrality measurements.

However, in most circumstance, the function is called by \link{cepa.all}.
}
\value{
    A \code{cepa.all} class object
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}
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
res.gsa = cepa.univariate.all(mat = eset, label = label, pc = PID.db$NCI)
}
}