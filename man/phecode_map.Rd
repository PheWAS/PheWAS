\name{phecode_map}
\alias{phecode_map}
\docType{data}
\title{
Dataframe containing phecode mappings for ICD9CM and ICD10CM.
}
\description{
This data frame maps every ICD9CM and ICD10CM code to its direct phecode.
This mapping is one to one.
}
\usage{data(phemap)}
\format{
  A data frame with 84532 observations on the following 3 variables.
  \describe{
    \item{\code{vocabulary_id}}{Character vector representing the vocabulary of the code- either 'ICD9CM' or 'ICD10CM'}
    \item{\code{code}}{Character vector representing the specific code to be mapped}
    \item{\code{phecode}}{Character vector representing the direct phecode mapping}
  }
}

\keyword{datasets}

\seealso{
\code{\link[PheWAS:phecode_rollup_map]{PheWAS::phecode_rollup_map}}
\code{\link[PheWAS:mapCodesToPhecodes]{mapCodesToPhecodes}}
}