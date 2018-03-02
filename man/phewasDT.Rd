\name{phewasDT}
\alias{phewasDT}

\title{
Wrapper to create a browsable datatable
}
\description{
This function returns a DT from the pacakge DT that has commonly used columns in an easy to browse format.
}
\usage{
phewasDT(results, addPhecodeInfo=T, OR.digits=2, p.digits=3, ...)
}
\arguments{
  \item{results}{
The results data frame from \code{phewas}.
}
  \item{addPhecodeInfo}{
Add the phecode information? Default is TRUE, which works with a phecode PheWAS. FALSE requires the inclusion of a "description" column.
}
  \item{OR.digits}{
Digits after the decimal place for the odds ratio (passed as the digits parameter to \code{round}). Default is 2.
}
  \item{p.digits}{
Digits after the decimal place for the p-value (passed as the digits parameter to \code{round}). Default is 3.
}
  \item{\dots}{
  Parameters to pass to \code{\link[DT:datatable]{datatable}}
}
}
\details{
This function is a short wrapper to create an easily browsable datatable using the package DT. It has few options, but this function serves as an instructive example as the interior calls of the function can be easily adapted to suit other purposes.
}

\author{
Robert Carroll
}

\seealso{
\code{\link[PheWAS:createPhewasTable]{createPhewasTable}} 
}
\examples{
\donttest{
#Generate some example data
ex=generateExample(hit="335")
#Extract the two parts from the returned list
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
#Create the PheWAS code table- translates the icd9s, adds exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count)
#Run the PheWAS
results=phewas(phenotypes,genotypes,cores=4)
#Create a quick datatable for results review
phewasDT(results)
}
}
\keyword{ utilities }

