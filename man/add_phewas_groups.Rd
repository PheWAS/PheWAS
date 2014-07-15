\name{add_phewas_groups}
\alias{add_phewas_groups}
\title{
Add PheWAS code groups to existing data.
}
\description{
\code{add_phewas_groups} adds the PheWAS code groups, found in the data \code{code_group_ranges}. Requires the full three digit codes in character format. This function is primarily to assist in plotting PheWAS results.
}
\usage{
add_phewas_groups(d, keep_newlines = T)
}
\arguments{
  \item{d}{
Data frame containing the attribute "phenotype" which has PheWAS codes.
}
  \item{keep_newlines}{
There are newlines in some longer group names to help with plotting which are kept by default. False removes these newlines.
}
}
\value{
Data frame with added columns: 
\item{groupnum}{The number of the group, ordered by PheWAS codes.}
\item{group}{The name of the group.}
\item{color}{The default color for the group when plotting.}
}

\keyword{ utilities }
