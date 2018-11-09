\name{phecode_rollup_map}
\alias{phecode_rollup_map}
\docType{data}
\title{
Dataframe containing phecode rollup mappings.
}
\description{
This data frame maps every phecode to itself and all its parent phecodes.
This mapping is one to many.
}
\usage{data(phemap)}
\format{
  A data frame with 3227 observations on the following 2 variables.
  \describe{
    \item{\code{code}}{a character vector representing the phecode to be rolled up}
    \item{\code{phecode_unrolled}}{a character vector representing the phecode and each of its parent codes}
  }
}

\keyword{datasets}

\seealso{
\code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}}
\code{\link[PheWAS:mapCodesToPhecodes]{mapCodesToPhecodes}}
}