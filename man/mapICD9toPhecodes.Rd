\name{mapICD9ToPhecodes}
\alias{mapICD9ToPhecodes}
\title{
Map ICD9 codes to phecodes
}
\description{
This function takes a vector of International Classification of Disease 9th Edition Clinical Modification codes and returns the phecodes represented. This is a many to many mapping.
}
\usage{
mapICD9ToPhecodes(..., icd9s, add.info=F, keep.icd9=F)
}
\arguments{
  \item{...}{
a data frame or arguments of either the form value or tag = value.
}
  \item{icd9s}{
a character vector of ICD9 billing codes, including trailing '0's.
}
  \item{add.info}{
a boolean informing if the returned data frame should include phecode descriptions and groups.
}
  \item{keep.icd9}{
a boolean requesting to keep the original ICD9 code column in the returned data frame.
}
}
\details{
By supplying ids and counts, one can convert data for multiple records at once, preserving count information. Note that this function does not collapse duplicate mappings; id-phecode pairs will not necessarily be unique.
}
\value{
A data frame containing the input columns in ... and the phecodes. The origial ICD9s and description are included if requested.
\item{...}{Other rows supplied e.g., ids or counts}
\item{icd9 }{Input ICD9}
\item{phecode }{phecodes}
\item{phecode_description}{The description of the phecode.}
\item{phecode_group}{The name of the phecode group.}
\item{phecode_groupnum}{The assigned number of the phecode group.}
}
\author{
Robert Carroll
Laura Wiley
}

\examples{
icd9codes=c("008","714.1","250.00")
phecodes=mapICD9ToPhecodes(icd9s=icd9codes)$phecode
head(phecodes)
}
\keyword{ utilities }
