\name{restrictPhewasByGender}
\alias{restrictPhewasByGender}
\title{
Add PheWAS code descriptions to existing data.
}
\description{
\code{restrictPhewasByGender} alters a PheWAS table, as from \code{\link[PheWAS:createPhewasTable]{createPhewasTable}}, to exclude individuals with non-applicable genders from certain phenotypes.
}
\usage{
restrictPhewasByGender(phenotypes,id.gender)
}
\arguments{
  \item{phenotypes}{
  The PheWAS table to have restrictions applied. The first column should be the id.
}
  \item{id.gender}{
  A data frame with the first column being the id and the second the gender, "M" or "F", of the individual. Individuals with any other specification will have all gender specific phenotypes set to NA.
}
}
\value{
The \code{phenotypes} data frame with NA values for individuals that do not match the gender for gender-specific codes.
}

\keyword{ utilities }
