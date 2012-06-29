\name{set.pathway.catalogue}
\alias{set.pathway.catalogue}
\title{store pathway data and pre-processing}
\description{
store pathway data and pre-processing
}
\usage{
set.pathway.catalogue(pathList, interactionList, mapping,
    min.node = 5, max.node = 500, min.gene = min.node, max.gene = max.node,
    ...)    
}
\arguments{
  \item{pathList}{list of pathways}
  \item{interactionList}{list of interactions}
  \item{mapping}{a data frame or matrix providing mappings from gene id to pathway node id. 
  The first column is node id and the second column is gene id.}
  \item{min.node}{minimum number of nodes in each pathway}
  \item{max.node}{maximum number of nodes in each pathway}
  \item{min.gene}{minimum number of genes in each pathway}
  \item{max.gene}{maximum number of genes in each pathway}
  \item{...}{other arguments}
}
\details{
Basicly, a pathway contains a list of interactions. The pathList argument
is a list where elements in the list is the vector of interaction IDs 
in the pathway. The interactions in the pathway can be got from a interaction
list pool represented as interactionList argument. The interactionList argument
stores the total interaction list in the pathway catalogue. It represents
as a three columns data frame or matrix where the first column is the interaction id,
the second column is the input node id and the third column is the output node id.

The mapping data frame provide the mapping from node id to gene id. The 
first column is the node id and the second column is the gene id.

}
\value{
    A \code{pathway.catalogue} class object
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{
data(PID.db)

catalogue = set.pathway.catalogue(pathList = PID.db$NCI$pathList[1:20],
            interactionList = PID.db$NCI$intertionList, mapping = PID.db$NCI$mapping)
}
}
