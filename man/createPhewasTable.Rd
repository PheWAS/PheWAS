\name{createPhewasTable}
\alias{createPhewasTable}
\title{
Creates a phenotype table from id, ICD-9 (or phewas) code, count data.
}
\description{
This function takes a data frame with three columns: id, icd9 codes, and counts. It returns a wide table with ICD-9 codes translated into phewas codes. It can optionally use the PheWAS exclusion criteria. Results can be formatted as true/false for a specified threshold or as summed continuous counts.
}
\usage{
createPhewasTable(id.icd9.count, min.code.count = 2, add.exclusions = T,
  translate = T, aggregate.fun = sum, id.gender)
}

\arguments{
  \item{id.icd9.count}{
Data frame with three columns: id, ICD-9 character string, and count.
}
  \item{min.code.count}{
The minimum code count to be considered a case. NA results in a continuous output.
}
  \item{add.exclusions}{
Apply PheWAS exclusions to the phewas codes.
}
  \item{translate}{
Should the input ICD-9 codes be translated to phewas codes?
}
  \item{aggregate.fun}{
Aggregate function for duplicated ICD9 or phewas codes in an individual.
}
  \item{id.gender}{
If supplied, restrict the phenotypes by gender. This should be a data frame with the first column being the id and the second the gender, "M" or "F", of the individual. Individuals with any other specification will have all gender specific phenotypes set to NA.
}
}
\details{
By default, this function returns a wide format data frame with boolean phenotypes suitable for PheWAS analysis. Specifying a \code{min.code.count=NA} will permit continuous code count phenotypes.

The default exclusions can be skipped with \code{add.exclusions=F}. In conjuntion with \code{translate=F} (and optionally adjusting \code{min.code.count} and \code{aggregate.fun}), one can use this function as a simple reshaping wrapper.
}
\value{
A data frame. The first column, named "id", contains the supplied id for each individual. The following columns are all present phewas codes. They contain T/F/NA for case/control/exclude or continuous/NA if min.code.count was NA.
}
\author{
Robert Carroll
Laura Wiley
}
\examples{
#Simple example
id_icd9_count=data.frame(id=c(1,2,2,2,3),icd9=c("714","250.11","714.1","714","250.11"),
  count=c(1,5,1,1,0))
createPhewasTable(id_icd9_count)
\donttest{
#Complex example
ex=generateExample(n=500,hit="335")
#Extract the icd9 data
id.icd9.count=ex$id.icd9.count
#Create the PheWAS code table
phenotypes=createPhewasTable(id.icd9.count)
phenotypes[1:10,1:10]
}
}
\keyword{ utilities }
