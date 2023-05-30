\name{restrictPhecodesBySex}
\alias{restrictPhecodesBySex}
\alias{restrictPhecodesByGender}
\title{
Restrict individuals for sex-specific phecodes 
}
\description{
\code{restrictPhecodesBySex} alters a table for PheWAS with phecodes, as from \code{\link[PheWAS:createPhenotypes]{createPhenotypes}}, to exclude individuals with non-applicable sexes from certain phenotypes.
}
\usage{
restrictPhecodesBySex(phenotypes,id.sex,sex.restriction=PheWAS::sex_restriction)
}
\arguments{
  \item{phenotypes}{
  The PheWAS table to have restrictions applied. The first column should be the id.
}
  \item{id.sex}{
  A data frame with the first column being the id and the second the sex, "M" or "F", of the individual. Individuals with any other specification will have all sex-specific phenotypes set to NA.
}
  \item{sex.restriction}{
  List that specifies sex-specific phecodes. By default uses the provided PheWAS::sex_restriction.
}
}
\value{
The \code{phenotypes} data frame with NA values for individuals that do not match the sex for sex-specific codes.
}

\keyword{ utilities }
