\name{cepa}
\alias{cepa}
\title{Access the significance of a pathway}
\description{
Using network centrality as weight to access the significance of a pathway.
Pathway nodes are taken as network units instead of genes.
}
\usage{
cepa(dif, bk, pathway, mapping = cbind(bk, bk), cen = "equal.weight",
     cen.name = if(is.function(cen)) deparse(substitute(cen))
                else if(mode(cen) == "name") deparse(cen)
                else cen,
     iter = 1000)        
}
\arguments{
  \item{dif}{differential gene list}
  \item{bk}{background gene list}
  \item{pathway}{\link[igraph]{igraph} object or edge list}
  \item{mapping}{a data frame or matrix providing mappings from gene id to pathway node id. 
  The first column is node id and the second column is gene id.}
  \item{cen}{centrality measuments, it can ce a string, function, 
  or function that has been quoted}
  \item{cen.name}{centrality measurement names. This argument should be set
  if the \code{cen} is a function. }
  \item{iter}{number of simulations}
}
\details{
The differential gene list and the background gene list should be indicated
with the same identifiers (e.g. gene symbol or refseq ID). All genes in
the differential gene list should exist in the background gene list.

The pathway should be expressed as a edge list represented as a matrix
or a \link[igraph]{igraph} object which has already embed the edge list.
Vertex in the pathway network should be indicated with node id. If every
gene in the pathway is a single node, then the node id is just the gene id.

The mapping data frame provide the mapping from node id to gene id. The 
first column is the node id and the second column is the gene id. By default,
it assumes every gene is the pathway is a single node.

If the centrality measurement is set as a string, only pre-defined "equal.weight",
"in.degree", "out.degree", "degree", "betweenness", "in.reach", "out.reach",
"reach", "in.spread", "out.spread" and "spread" are allowed. More centrality
measurements can be used by setting it as a function (such as closeness,
cluster coefficient).

}
\value{
A \code{cepa} class object which is a list with 12 components:
  \item{score}{Pathway score}
  \item{ds}{four descriptive values of centralities in the pathway
  (max, 75th quantile, median, min)}
  \item{ds.simulation}{four descriptive values of centralities in simulations}
  \item{p.value}{p-value calculated by CePa}
  \item{p.ora}{p-value calculated by ORA}
  \item{simulation}{pathway scores in simulations}
  \item{centrality}{centrality measurement name}
  \item{weight}{weight of nodes in the pathway}
  \item{is.dif.node}{logical vector identifying whether the node is differentially affected}
  \item{node.name}{names of nodes in pathway. Node name is composed by the names of the member genes}
  \item{count}{four count of genes and nodes in the pathway}
  \item{pathway}{pathway in \link[igraph]{igraph} class}
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
######################################
# a simulation
######################################
require(igraph)
path.genes = 0:199
bk = 0:9999
pathway = barabasi.game(200)
# in the pathway, there are 40 differential genes
dif = c(sample(0:199, 40), 200:1159)

res = cepa(dif, bk, pathway, cen="degree", iter = 500)
res

#####################################
# a real pathway
#####################################
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
