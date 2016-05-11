\name{addPhecodeInfo}
\alias{addPhecodeInfo}
\alias{addPhewasDescription}
\alias{addPhewasGroups}
\title{
Add phecode descriptions and group information to existing data.
}
\description{
\code{addPhewasDescription} adds the phecode description and group information, found in the data frame \code{pheinfo}. Requires the full three digit phecodes in character format. This function replaces addPhewasDescription and addPhewasGroups.
}
\usage{
addPhecodeInfo(data, descriptions=T, groups=T, groupnums=F)
}
\arguments{
  \item{data}{
Data frame containing the a column with phecodes. Can additionally be a character vector of phecodes. See details for requirements.
}
  \item{descriptions}{
Add the phecode descriptions? Default is TRUE (yes).
}
  \item{groups}{
Add the phecode group names? Default is TRUE (yes).
}
  \item{groupnums}{
Add the phecode group numbers? Default is FALSE (no). Used in plotting
}
}
\details{
This function provides a quick way to add phecode descriptions and group names.
It will select the first column matching any of "phecode", "phewas" or "pheno" to merge on. If none is identified, it reports a warning and attempts the merge with the first column.
Phecodes must be character vectors, otherwise they may have inaccurate mappings. If PheWAS codes are factors it will convert them to characters, give a warning, and attempt to map them.
If a character (or factor) variable or vector is passed in, it will return code descriptions for those codes.
}
\value{
Data frame with added columns: 
\item{phecode_description}{The description of the phecode.}
\item{phecode_group}{The name of the phecode group.}
\item{phecode_groupnum}{The assigned number of the phecode group.}
}

\keyword{ utilities }
