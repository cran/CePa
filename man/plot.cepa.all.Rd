\name{plot.cepa.all}
\alias{plot.cepa.all}
\title{plot the cepa.all object}
\description{
plot the cepa.all object
}
\usage{
\method{plot}{cepa.all}(x, id = NULL, cen = 1, type = c("graph", "null"),
    node.name = NULL, node.type = NULL,
    adj.method = "none", only.sig = FALSE,
    cutoff = ifelse(adj.method == "none", 0.01, 0.05), ...)                    
}
\arguments{
  \item{x}{a \link{cepa.all} object}
  \item{id}{index or the name for the pathway}
  \item{cen}{index or the name for the centrality}
  \item{type}{If the aim is to plot single pathway, then this argument is to identify the kind of the plotting.}
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
x, adj.method, only.sig, cutoff arguments. 

Second, it can draw figures of single
pathway under specific centrality measurement. Under this circumstance, 
this function is just a wrapper of \link{plot.cepa}. To do it, 
users should set x, id, cen, type, node.name and node.type arguments. It must
be noted that these two kinds of arguments should not be mixed.

}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{
data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
plot(res.ora)
plot(res.ora, id = 3)
plot(res.ora, id = 3, type = "null")

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
plot(res.gsa)
plot(res.gsa, id = 3, cen = 2)
plot(res.gsa, id = 3, cen = 2, type = "null")
}
}
