\name{mapPhecodesToExclusions}
\alias{mapPhecodesToExclusions}
\title{
Map phecodes to their exclusions
}
\description{
This function maps phecodes (optionally with ids for individuals) to a set of phecode exclusions.
}
\usage{
mapPhecodesToExclusions(phecodes, ids)
}
\arguments{
  \item{phecodes}{
A vector of phecodes.
}
  \item{ids}{
An optional vector of ids to pair with the provided phecodes.
}
}

\value{
A data frame containing phecodes and their exclusions. IDs for those codes and exclusions are included if they were supplied.
\item{id}{If ids were provided, the individual ids are included as the first column}
\item{exclusion_criteria}{Input phecodes}
\item{exclusion}{The exclusion phecodes for the codes provided}
}
\author{
Robert Carroll
Laura Wiley
}

\keyword{ utilities }
