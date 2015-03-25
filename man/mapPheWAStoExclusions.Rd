\name{mapPheWAStoExclusions}
\alias{mapPheWAStoExclusions}
\title{
Map PheWAS codes to their exclusions
}
\description{
This function maps phewas codes (optionally with ids for individuals) to a set of PheWAS code exclusions.
}
\usage{
mapPheWAStoExclusions(phewas.codes, ids)
}
\arguments{
  \item{phewas.codes}{
A vector of PheWAS codes.
}
  \item{ids}{
An optional vector of ids to pair with the provided PheWAS codes.
}
}

\value{
A data frame containing phewas codes and their exclusions. IDs for those codes and exclusions are included if they were supplied.
\item{id}{If ids were provided, the individual ids are included as the first column}
\item{exclusion_criteria}{Input PheWAS codes}
\item{exclusion}{The exclusion PheWAS codes for the codes provided}
}
\author{
Robert Carroll
Laura Wiley
}

\keyword{ utilities }
