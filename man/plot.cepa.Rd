\name{plot.cepa}
\alias{plot.cepa}
\title{Plot the cepa object}
\description{
Plot the cepa object
}
\usage{
\method{plot}{cepa}(x, node.name = NULL, node.type = NULL,
    graph.node.max.size = 20, graph.node.min.size = 3, graph.layout.method = NULL, ...)                 
}
\arguments{
  \item{x}{\link{cepa} object}
  \item{node.name}{node.name for each node}
  \item{node.type}{node.type for each node}
  \item{graph.node.max.size}{max size of the node in the graph}
  \item{graph.node.min.size}{min size of the node in the graph}
  \item{graph.layout.method}{function of the layout method. For the list
  of available methods, see \link[igraph]{layout}}
  \item{...}{other arguments}
}
\details{
There are four figures in the plotting.
A) Distribution of reach centrality in the pathway under simulation. 
The distribution of node centralities in each simulation is measures by maximum value, 
the 75th quantile, median value and minimum value; B) Distribution of centrality 
of all nodes in the real pathway; C) Histogram of simulated scores in the pathway; 
D) Graph view of the pathway where the size of node is proportional to centrality 
value of the node. In figure A and B, dot on the plot is randomly shifted to make
a more clear visualization of the distribution of data. In figure A and C, 
the real pathway is marked with red line.
 
The function returns a \link[igraph]{igraph} object, so if users don't satisfy 
with the network visualization, they can draw graphs of the network with their
own settings.

}
\value{
A \link[igraph]{igraph} object of the pathway
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
plot(res)
}
