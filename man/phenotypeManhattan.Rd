\name{phenotypeManhattan}
\alias{phenotypeManhattan}
\title{
Plot phenotype association data in a Manhattan plot.
}
\description{
\code{phenotype_manhattan} wraps \code{phenotypePlot} to provide some intuitive Manhattan plot centric options. Can plot data directly from \code{phewas}.
}
\usage{
phenotypeManhattan(d, suggestive.line = 0.05, genomewide.line, sort.by.p = F, 
  sort.by.category.p = F, OR.size = F, OR.direction = F, annotate.level, 
  y.axis.interval = 5, y.axis.label = expression(-log[10](italic(p))), max.y, ...)
}
\arguments{
  \item{d}{
A data frame containing \code{phenotype} and \code{p}.
}
  \item{suggestive.line}{
Draw a blue line at this p-value to show a suggestive significance threshold. Defaults to 0.05. An NA value will result in no line.
}
  \item{genomewide.line}{
Draw a red line at this p-value to show a genome wide significance threshold. Defaults to suggestive.line divided by the number of non-NA p-values. An NA value will result in no line.
}
  \item{sort.by.p}{
Sort the plot by the p-value. Defaults to false.
}
  \item{sort.by.category.p}{
Sort the plot by pvalue, within each phenotype group. Requires group attributes in \code{d}. Default is FALSE.
}
  \item{OR.size}{
Adjust point size based on odds ratios? Requires d to contain the \code{OR} column. Default is FALSE.
}
  \item{OR.direction}{
Adjust point shape based on odds ratio direction? Requires d to contain the \code{OR} column. Default is FALSE.
}
  \item{annotate.level}{
Annotate points at or below this \code{p}. Defaults to genomewide.line.
}
  \item{y.axis.interval}{
The -log10 interval for y axis labeling. Defaults to 5.
}
  \item{y.axis.label}{
The label for the y axis. Defaults to a stylized "-log10(p)".
}
  \item{max.y}{
The maximum \code{-log10(p)} on the axis.
}
  \item{\dots}{
Further parameters to be passed to \code{phenotype_plot}.
}
}
\details{
Produces a Manhattan plot by phenotype. Further options are available from the inner function, \code{phenotypePlot}.
}
\value{
A \code{ggplot2} plot with phenotypes against their negative log transformed p-values.
}

\note{
These parameters map to ones in \code{phenotype_plot}, with different names (and some extra processing) based on context (i.e., generic value plotting versus association study plotting).
}

\seealso{
This function can be called from \code{\link[PheWAS:phewasManhattan]{phewasManhattan}} if all phenotypes are PheWAS codes.
For further parameters, see \code{\link[PheWAS:phenotypePlot]{phenotypePlot}}.
}


\keyword{ hplot }
