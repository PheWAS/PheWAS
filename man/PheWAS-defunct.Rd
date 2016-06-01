\name{PheWAS-defunct}
\alias{PheWAS-defunct}
\title{
Defunct methods in the PheWAS package
}
\description{
List of methods that have been removed from the PheWAS package over time.
}

\details{
Several changes have been made to promote consistency with the newer language around phecodes: 
\itemize{
\item \code{addPhewasDescription} and 
\item \code{addPhewasGroups} have been replaced by \code{\link[PheWAS:addPhecodeInfo]{addPhecodeInfo}}. 
\item \code{mapICD9stoPheWAS} has been replaced with \code{\link[PheWAS:mapICD9toPhecodes]{mapICD9toPhecodes}}.
\item \code{mapPheWAStoExclusions} has been replaced with \code{\link[PheWAS:mapPhecodesToExclusions]{mapPhecodesToExclusions}}.
}

}

\seealso{
\code{\link[PheWAS:addPhecodeInfo]{addPhecodeInfo}}
\code{\link[PheWAS:mapICD9toPhecodes]{mapICD9toPhecodes}}
\code{\link[PheWAS:mapPhecodesToExclusions]{mapPhecodesToExclusions}}
}
