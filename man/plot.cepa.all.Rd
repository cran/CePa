\name{plot.cepa.all}
\alias{plot.cepa.all}
\title{plot the cepa.all object}
\description{
plot the cepa.all object
}
\usage{
\method{plot}{cepa.all}(x, id = NULL, cen = 1,
    node.name = NULL, node.type = NULL,
    adj.method = "none", only.sig = FALSE,
    cutoff = ifelse(adj.method == "none", 0.01, 0.05), ...)                    
}
\arguments{
  \item{x}{\link{cepa.all} object}
  \item{id}{index or the name for the pathway}
  \item{cen}{index or the name for the centrality}
  \item{node.name}{node.name for each node}
  \item{node.type}{node.type for each node}
  \item{adj.method}{method of \link[stats]{p.adjust}}
  \item{only.sig}{whether show all pathways}
  \item{cutoff}{cutoff for significance}
  \item{...}{other arguments}
}
\details{
This function has two applications. First, it can draw heatmaps of p-values
of all pathways under different centrality measurements. To do it, users should set
x, adj.method, only.sig, cutoff arguments. Second, it can draw figures of single
pathway under specific centrality measurement. Under this circumstance, 
this function is just a wrapper of \link{plot.cepa}. To do it, 
users should set x, id, cen, node.name and node.type arguments. It must
be noted that these two kinds of arguments should not be mixed.

}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{
data(pathway.db)
data(gene.list)
res = cepa.all(gene.list$dif, gene.list$bk, pathway.db$NCI$pathList,
               pathway.db$NCI$interactionList, pathway.db$NCI$mapping)
plot(res, adj.method = "BH")
plot(res, id = 1, cen = "in.degree")
}
}
