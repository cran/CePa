\name{report}
\alias{report}
\title{Generate report for CePa analysis}
\description{
Generate report for CePa analysis
}
\usage{
report(x, adj.method = "none", cutoff = ifelse(adj.method == "none", 0.01, 0.05), 
       template.file = paste(system.file(package = "CePa"), "/extdata/cepa.template", sep=""),
       only.sig = TRUE, ...)
}
\arguments{
  \item{x}{\link{cepa.all} object}
  \item{adj.method}{methods in \link[stats]{p.adjust}}
  \item{cutoff}{cutoff for significance}
  \item{template.file}{path of the template file}
  \item{only.sig}{whether to generate detailed report for every pathway.
  If it is set to FALSE, the page for every pathway under every centrality
  would be generated.}
  \item{...}{other arguments}
}
\details{
    The report would locate at the current working directory. View the report
	by clicking \code{index.html} in the report directory.
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa.all}}
}
\examples{
\dontrun{
data(PID.db)
data(gene.list)
res = cepa.all(gene.list$dif, gene.list$bk, PID.db$NCI$pathList,
               PID.db$NCI$interactionList, PID.db$NCI$mapping)
report(res)
}
}
