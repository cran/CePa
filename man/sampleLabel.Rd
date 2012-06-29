\name{sampleLabel}
\alias{sampleLabel}
\title{Generate data structure of sample labels}
\description{
Generate data structure of sample labels
}
\usage{
sampleLabel(label, treatment, control)
}
\arguments{
  \item{label}{sample label vector}
  \item{treatment}{treatment label}
  \item{control}{control label}
}
\value{
A \code{sampleLabel} class object
}
\author{Zuguang Gu}
\examples{
sampleLabel(c("A", "B", "B", "A", "A", "A", "B", "B"), treatment = "A", control = "B")
}
