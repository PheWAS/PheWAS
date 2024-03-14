% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restrictPhecodesBySex.R
\name{restrictPhecodesBySex}
\alias{restrictPhecodesBySex}
\title{Add PheWAS code descriptions to existing data.}
\usage{
restrictPhecodesBySex(phenotypes, id.sex)
}
\arguments{
\item{phenotypes}{The PheWAS table to have restrictions applied. The first
column should be the id.}

\item{id.sex}{A data frame with the first column being the id and the second
the gender, "M" or "F", of the individual. Individuals with any other

specification will have all gender specific phenotypes set to NA.}
}
\value{
The \code{phenotypes} data frame with NA values for individuals that
do not match the sex for sex-specific codes.
}
\description{
\code{restrictPhecodesBySex} alters a table for PheWAS with phecodes, as from
\code{\link[PheWAS:createPhenotypes]{createPhenotypes}}, to exclude
individuals with non-applicable sexes from certain phenotypes.
}
\examples{
data <- sample_data
restrictPhecodesBySex(data$id.vocab.code.count, data$id.sex)
}