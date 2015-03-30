\name{phenotypePlot}
\alias{phenotypePlot}

\title{
Plot phenotypes with annotation.
}
\description{
\code{phenotypePlot} plots simple phenotype data, useful in circumstances where p-values are not the measure under study, as in \code{\link[PheWAS:phewasManhattan]{phewasManhattan}}.
}
\usage{
phenotypePlot(d, max.y, max.x, suggestive.line, significant.line,
  size.x.labels = 9, size.y.labels = 9, switch.axis = F, sort.by.value = F,
  sort.by.category.value = F, annotate.phenotype.description, 
  annotate.angle = 30, annotate.size = 5, annotate.level, 
  annotate.phenotype = F, annotate.snp.w.phenotype = F, annotate.snp = F,
  annotate.snp.angle = 0, annotate.list, annotate.only.largest = T,
  lc.labels = F, x.group.labels = T, x.phenotype.labels = F, sizes = F,
  direction = F, point.size = 3, use.color = T, color.palette, 
  title = paste0("Phenotype Plot ", date()), x.axis.label = "Phenotypes", 
  y.axis.label = "Values", y.axis.interval = 10)
}
\arguments{
  \item{d}{
A data frame containing at least the columns "value" and "phenotype".
}
  \item{max.y}{
Maximum \code{value} for the plot. If missing, defaults to the maximum value.
}
  \item{max.x}{
Maximum number of \code{phenotype}s for the plot. If missing, defaults to the total number.
}
  \item{suggestive.line}{
Blue line drawn at \code{value} to emphasize interesting data points. Missing or an NA value will result in no line.
}
  \item{significant.line}{
Red line drawn at \code{value} to emphasize interesting data points. Missing or an NA value will result in no line.
}
  \item{size.x.labels}{
The size of the x axis labels.
}
  \item{size.y.labels}{
The size of the y axis labels.
}
  \item{switch.axis}{
Switch the X and Y axis? Currently non-functioning. Defaults to false.
}
  \item{sort.by.value}{
Sort the plot by highest to lowest \code{value}? Defaults to false.
}
  \item{sort.by.category.value}{
Sort the plot by \code{value} within each phenotype group. Requires group attributes in \code{d}. Defaults to false.
}
  \item{annotate.phenotype.description}{
Contains either TRUE or a data frame that contains description annotations for each phenotype. Should contain columns "phenotype" and "description". Missing or FALSE yields no description annotation.
}
  \item{annotate.angle}{
Angle for annotation text.
}
  \item{annotate.size}{
Size of annotation text
}
  \item{annotate.level}{
Minimum \code{value} to annotate points.
}
  \item{annotate.phenotype}{
Include the phenotype name in the annotation? Default is FALSE.
}
  \item{annotate.snp.w.phenotype}{
Include the snp in the phenotype annotation? Default is FALSE.
}
  \item{annotate.snp}{
Annotate the SNP? This is a second annotation, default is FALSE.
}
  \item{annotate.snp.angle}{
Angle of annotation for the SNP only annotation.
}
  \item{annotate.list}{
A vector of \code{phenotype}s to force annotation, regardless of significance.
}
  \item{annotate.only.largest}{
Should only the largest point for each phenotype be annotated? Default is TRUE, which only shows text for each phenotype once, even if multiple points pass the annotation threshold or are listed.
}
  \item{lc.labels}{
Force the labels to lower case? Default is FALSE.
}
  \item{x.group.labels}{
Label the phenotype groups? Requires group attributes in \code{d}. Default is TRUE.
}
  \item{x.phenotype.labels}{
Label every phenotype? Default is no, cannot be used in conjunction with \code{x.group.labels}.
}
  \item{sizes}{
Adjust point size based on \code{size} in \code{d}? Default is FALSE.
}
  \item{direction}{
Adjust point shape based on \code{direction} in \code{d}? Default is FALSE.
}
  \item{point.size}{
Size of the points. Default is 2.
}
  \item{use.color}{
Color the points? Requires \code{color.palette} or \code{color} in \code{d}. Default is TRUE. 
}
  \item{color.palette}{
An alternate color palette. Requires at least the number of groups that exist.
}
  \item{title}{
Title for the plot.
}
  \item{x.axis.label}{
Label for the \code{phenotype} axis.
}
  \item{y.axis.label}{
Label for the \code{value} axis.
}
  \item{y.axis.interval}{
Spacing for the ticks on the \code{value} Axis.
}
}
\value{
A ggplot2 plot object with the requested phenotype plot.
}
\references{
This code is a heavy modification and update from Stephen Turner's code at http://gettinggeneticsdone.blogspot.com/2010/01/gwas-manhattan-plots-and-qq-plots-using.html
}

\seealso{
\code{\link[PheWAS:phewasManhattan]{phewasManhattan}} and \code{\link[PheWAS:phenotypeManhattan]{phenotypeManhattan}} for wrappers that permit simpler plotting of common plots.
}

\keyword{ hplot }

