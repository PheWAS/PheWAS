\name{code_group_ranges}
\alias{code_group_ranges}
\docType{data}
\title{
Data frame containing phewas code group information.
}
\description{
This data frame contains the name, PheWAS code ranges, and default plotting color for the 14 PheWAS code groups.
}
\usage{data(code_group_ranges)}
\format{
  A data frame with 14 observations on the following 5 variables.
  \describe{
    \item{\code{groupnum}}{integer vector with the group number}
    \item{\code{min}}{a character vector with the smallest code in each group}
    \item{\code{max}}{a character vector with the largest possible code in each group}
    \item{\code{group}}{a character vector containing the name of the group, suitable for plotting}
    \item{\code{color}}{a character vector containing the default color for this group in phewas manhattan plots}
  }
}
\source{
http://phewascatalog.org
}
\keyword{datasets}
