\name{phewas_exclude}
\alias{phewas_exclude}
\docType{data}
\title{
Dataframe containing phewas code mappings to their exclusions
}
\description{
This data frame maps phewas codes to their exclusions: analyses in which an individual with a particular code may not be an appropriate control.
This mapping is many to many.
}
\usage{data(phewas_exclude)}
\format{
  A data frame with 30033 observations on the following 2 variables.
  \describe{
    \item{\code{code}}{A character vector of PheWAS codes.}
    \item{\code{exclusion_criteria}}{A character vector of PheWAS codes; an individual having this code would exclude one from an analysis of the PheWAS code \code{code}.}
  }
}
\source{
http://phewascatalog.org
}
\keyword{datasets}
