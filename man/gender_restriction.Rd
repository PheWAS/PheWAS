\name{gender_restriction}
\alias{gender_restriction}
\docType{data}
\title{
Data frame containing potential gender-restricted phewas codes.
}
\description{
Some phewas codes are only typical for one gender or the other; this data frame describes those restrictions.
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
