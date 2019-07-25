\name{phecode_map_icd10}
\alias{phecode_map_icd10}
\docType{data}
\title{
Dataframe containing phecode mappings for ICD10.
}
\description{
This data frame maps each ICD10 code to its directly mapped phecode(s).
Can be provided to \code{\link[PheWAS:createPhenotypes]{createPhenotypes}} in the \code{vocabulary.map} parameter to convert ICD10 codes for use in the PheWAS methods.
It is the map 1.2b1.
\code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}} contains mapping for ICD9CM and ICD10CM (ie, US use cases).
}
\usage{data(phemap)}
\format{
  A data frame with 9505 observations on the following 3 variables.
  \describe{
    \item{\code{vocabulary_id}}{Character vector representing the vocabulary of the code- always 'ICD10'}
    \item{\code{code}}{Character vector representing the specific code to be mapped}
    \item{\code{phecode}}{Character vector representing the direct phecode mapping}
  }
}

\keyword{datasets}

\seealso{
\code{\link[PheWAS:phecode_rollup_map]{PheWAS::phecode_rollup_map}}
\code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}}
\code{\link[PheWAS:mapCodesToPhecodes]{mapCodesToPhecodes}}
}