\name{addPhewasDescription}
\alias{addPhewasDescription}
\title{
Add PheWAS code descriptions to existing data.
}
\description{
\code{addPhewasDescription} adds the PheWAS code description, found in the data \code{pheinfo}. Requires the full three digit codes in character format.
}
\usage{
addPhewasDescription(data, keep.unmatched.rows=F,for.plots=F)
}
\arguments{
  \item{data}{
Data frame containing the attribute with PheWAS codes. Can additionally be a character vector of PheWAS codes. See details for requirements.
}
  \item{keep.unmatched.rows}{
Keep rows that have no matching PheWAS descriptions? This is passed into keep.y in \code{merge}.
}
  \item{for.plots}{
Renames the output to include "phenotype" and "description". Default is \code{FALSE}.
}
}
\details{
This function provides a quick way to add PheWAS code descriptions.
It will select the first column matching either "phewas" or "pheno" to merge on. If neither are identified, it reports a warning and attempts the merge with the first column.
PheWAS codes must be character vectors, otherwise they may have inaccurate mappings. If PheWAS codes are factors it will convert them to characters, give a warning, and attempt to map them.
If a character (or factor) variable or vector is passed in, it will return codes descriptions for those codes.
}
\value{
Data frame with added columns: 
\item{phewas_code}{The number of the group, ordered by PheWAS codes.}
\item{phewas_description}{The name of the group.}
Note these names may be different if \code{for.plots=TRUE}.
}

\keyword{ utilities }
