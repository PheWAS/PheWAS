\name{phewasManhattan}
\alias{phewasManhattan}
\title{
Plot PheWAS results in a Manhattan plot.
}
\description{
This is a wrapper function that provides out of box capability to plot \code{phewas} function results.
}
\usage{
phewasManhattan(d, add.phewas.descriptions = T, ...)
}
\arguments{
  \item{d}{
Data frame containing \code{phenotype} and \code{p}.
}
  \item{add.phewas.descriptions}{
Annotate with PheWAS code descriptions? Default is TRUE.
}
  \item{\dots}{
Further parameters to be passed to \code{phenotypeManhattan} or \code{phenotypePlot}.
}
}
\details{
Simply plots a PheWAS Manhattan plot. Further options are available from the inner functions, \code{phenotypeManhattan} and \code{phenotypePlot}.
}
\value{
A \code{ggplot2} plot with PheWAS codes against their -log10 transformed p-values.
}

\seealso{
\code{\link[PheWAS:phenotypePlot]{phenotypePlot}} and \code{\link[PheWAS:phenotypeManhattan]{phenotypeManhattan}} for other available options.
\code{\link[PheWAS:PheWAS-package]{The PheWAS package}} for a complete example.
}
\keyword{ hplot }

