\name{PID.db}
\alias{PID.db}
\docType{data}
\title{pathway catalogues from Pathway Interaction Database(PID)}
\description{
pathway catalogues from Pathway Interaction Database(PID)
}
\usage{data(PID.db)}
\details{
The pathway data is parsed from XML format file provided by PID FTP site.
The perl code for parsing pathway data can be found at author's website
(\url{http://mcube.nju.edu.cn/jwang/lab/soft/cepa/index.html}).

There are four pathway catalogues which are NCI_Nature-curated, BioCarta, KEGG
and Reactome.

The node can be proteins, complex, families and none-gene nodes in real biological pathways. PID 
provides mapping from Uniprot ID to node id. In this package, gene symbol 
is selected as the gene identifier ID. The mapping data from gene symbol 
to UniProt ID is firstly get as gene symbol to HGNC id from HGNC and 
HGNC id to UniProt ID from UniProt database. 
}
\value{
A list containing four component:
\item{NCI}{NCI_Nature-curated pathway catalogue}
\item{BioCarta}{BioCarta pathway catalogue}
\item{KEGG}{KEGG pathway catalogue}
\item{Reactome}{Reactome pathway catalogue}

Each pathway catalogue is a \code{pathway.catalogue} class object. Each pathway
catalogue can be used directly in \link{cepa.all} and \link{cepa}

}
\source{
\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/NCI-Nature_Curated.xml.gz}

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/BioCarta.xml.gz} 

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/Reactome.xml.gz} 

\url{ftp://ftp1.nci.nih.gov/pub/PID/XML/KEGG.xml.gz}
}
\examples{
data(PID.db)
names(PID.db)
PID.db$NCI
plot(PID.db$NCI)
}