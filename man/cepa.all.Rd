\name{cepa.all}
\alias{cepa.all}
\title{Apply CePa algorithm on a list of pathways}
\description{
Apply CePa algorithm on a list of pathways
}
\usage{
cepa.all(dif, bk, pathList, interactionList, mapping = cbind(bk, bk),
         cen = c("equal.weight", "in.degree", "out.degree", "betweenness", "in.reach", "out.reach"),
         cen.name = sapply(cen, function(x){if(mode(x) == "name") deparse(x) else x}), iter = 1000,
         min.node = 5, max.node = 500)    
}
\arguments{
  \item{dif}{differential gene list}
  \item{bk}{background gene list}
  \item{pathList}{list of pathways}
  \item{interactionList}{list of interactions}
  \item{mapping}{a data frame or matrix providing mappings from gene id to pathway node id. 
  The first column is node id and the second column is gene id.}
  \item{cen}{centrality measuments, it can ce a string, or function has been quote}
  \item{cen.name}{centrality measurement names}
  \item{iter}{number of simulations}
  \item{min.node}{minimum number of nodes}
  \item{max.node}{maximum number of nodes}
}
\details{
The differential gene list and the background gene list should be indicated
with the same identifiers (e.g. gene symbol or refseq ID). All genes in
the differential gene list should exist in the background gene list.

Basicly, a pathway contains a list of interactions. The pathList argument
is a list where elements in the list is the vector of interaction IDs 
in the pathway. The interactions in the pathway can be got from a interaction
list pool represented as interactionList argument. The interactionList argument
stores the total interaction list in the pathway catalogue. It represents
as a three columns data frame or matrix where the first column is the interaction id,
the second column is the input node id and the third column is the output node id.

The mapping data frame provide the mapping from node id to gene id. The 
first column is the node id and the second column is the gene id. By default,
it assumes every gene is the pathway is a single node.

If the centrality measurement is set as a string, only pre-defined "equal.weight",
"in.degree", "out.degree", "degree", "betweenness", "in.reach", "out.reach",
"reach", "in.spread", "out.spread" and "spread" are allowed. More centrality
measurements can be used by setting it as a function (such as closeness,
cluster coefficient). In the cepa.all function, we recommand users choose
at least two centrality measurements.
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
data(gene.list)
res = cepa.all(gene.list$dif, gene.list$bk, PID.db$NCI$pathList,
               PID.db$NCI$interactionList, PID.db$NCI$mapping)
}
}
