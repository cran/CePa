\name{CePa-package}
\alias{CePa}
\alias{CePa-package}
\docType{package}
\title{Centrality-based pathway enrichment}
\description{
Use pathway topology information to assign weight to pathway nodes.}
\details{
Gene set enrichment analysis is broadly used in microarray data analysis
aimed to find which biological functions are affected by a group of 
related genes behind the massive information. A lot of methods have been
developed under the framework of over-represented analysis (ORA) such
as \code{GOstats} and \code{GSEABase}. For a specific
form of gene sets, biological pathways are collections of correlated genes/proteins, 
RNAs and compounds that work together to regulate specific biological
processes. Instead of just being a list of genes, a pathway contains 
the most important information that is how the member genes interact 
with each other. Thus network structure information is necessary for
the intepretation of the importance of the pathways.

In this package, the original pathway enrichment method
(ORA) is extended by introducing network centralities as the weight 
of nodes which have been mapped from differentially expressed genes 
in pathways. There are two advantages compared to former work.
First, for the diversity of genes' characters and the difficulties of 
covering the importance of genes from all aspects, we do not design a 
fixed measurement for each gene but set it as an optional parameter in the model. 
Researchers can select from candidate choices where different measurement 
reflects different aspect of the importance of genes. 
In our model, network centralities are used to measure the importance of genes in pathways. 
Different centrality measurements assign the importance to nodes from different aspects. 
For example, degree centrality measures the amount of neighbours that 
a node directly connects to, and betweenness centrality measures how many 
information streams must pass through a certain node. Generally speaking, 
nodes having large centrality values are central nodes in the network. 
It's observed that nodes represented as metabolites, proteins or genes 
with high centralities are essential to keep the steady state of biological networks. 
Moreover, different centrality measurements may relate to different biological functions. 
The selection of centralities for researchers depends on what kind of genes 
they think important. Second, we use nodes as the basic units of pathways 
instead of genes. We observe that nodes in the pathways include different 
types of molecules, such as single gene, complex and protein families. 
Assuming a complex or family contains ten differentially expressed member genes, 
in traditional ORA, these ten genes behave as the same position as other
genes represented as single nodes, and thus they have effect of ten. 
It is not proper because these ten genes stay in a same node in the 
pathway and make functions with the effect of one node. Also, 
a same gene may locate in different complexes in a pathway and if 
taking the gene with effect of one, it would greatly decrease the importance 
of the gene. Therefore a mapping procedure from genes to pathway nodes 
is applied in our model. What's more, the nodes in pathways also include 
none-gene nodes such as microRNAs and compounds. These nodes also 
contribute to the topology of the pathway. So, when analyzing pathways, 
all types of nodes are retained.

}
\author{
Zuguang Gu <jokergoo@gmail.com>
}
\keyword{package}
\examples{
\dontrun{
# the use of CePa package is simple
# the following is the standard process
# 1. load the pathway database
data(PID.db)

# 2. get your differential gene list
#    if you are using the PID.db pathway database
#    please ensure your gene id is official symbol (this is important)
# 3. also, your background gene list, should be as symbol
# 4. run CePa algorithm
#    here, dif.list is your differential gene list
#    and bk.list is your background gene list
res = cepa.all(dif.list, bk.list, PID.db$NCI$pathList,
               PID.db$NCI$interactionList, PID.db$NCI$mapping)
               
# 5. obtain the report of your analysis
report(res, adj.method = "BH", cutoff = 0.05, only.sig = FALSE)
# just wait for a moment, since there are a lot of figures to generate!
}
}
