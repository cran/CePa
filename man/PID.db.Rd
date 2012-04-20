\name{PID.db}
\alias{PID.db}
\docType{data}
\title{pathway catalogues from PID}
\description{
pathway catalogues from PID
}
\usage{data(PID.db)}
\details{
The pathway data is parsed from XML format file provided by PID FTP site.
The perl code for parsing pathway data can be found at author's website
(\url{http://mcube.nju.edu.cn/jokergoo/softwares.html}).

There are four pathway catalogues which are NCI_Nature-curated, BioCarta, KEGG
and Reactome.

The node can be proteins, complex, families and none-gene nodes. PID 
provides mapping from Uniprot ID to node id. In this package, gene symbol 
is selected as the gene identifies ID. The mapping data from gene symbol 
to UniProt ID is firstly get as gene symbol to HGNC id from HGNC and 
HGNC id to UniProt ID from UniProt database. 
}
\value{
A list containing four component:
\item{NCI}{NCI_Nature-curated pathway catalogue}
\item{BioCarta}{BioCarta pathway catalogue}
\item{KEGG}{KEGG pathway catalogue}
\item{Reactome}{Reactome pathway catalogue}

Foreach pathway catalogue, there are five components:

\item{pathList}{List, each item in the list is the vector of interaction ids}
\item{interactionList}{Data frame, the first column is the interaction id,
the second column is the input node id, the third column is the output
node id}
\item{mapping}{Data frame, the first column is the node id, the second
column is the gene id}
\item{node.name}{vector, names for the vector are the node ids}
\item{node.type}{vector, names for the vector are the node ids. Node types
contain protein, complex, compound and rna}
\item{created}{version of PID pathway data}
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
}