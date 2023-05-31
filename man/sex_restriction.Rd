\name{sex_restriction}
\alias{sex_restriction}
\docType{data}
\title{
Data frame containing potential sex-restricted phecodes.
}
\description{
Some phecodes are only typical for one sex or the other; this data frame describes those restrictions.
}
\usage{data(pheinfo)}
\format{
  A data frame with 1817 observations on the following 3 variables.
  \describe{
    \item{\code{phecode}}{The PheWAS code}
    \item{\code{male_only}}{Logical value for a male-only analysis}
    \item{\code{female_only}}{Logical value for a male-only analysis}
  }
}
\source{
http://phewascatalog.org
}
\keyword{datasets}
