\name{PheWAS_Plotting}
\alias{PheWAS_Plotting}
\alias{phewasManhattan}
\alias{phenotypeManhattan}
\alias{phenotypePlot}
\title{
Plotting methods for PheWAS results, phenotype association studies, or generic phenotype values
}
\description{
These functions permit the plotting of a variety of data types found when working with phenotypes. The functions are sequentially nested calls as outlined below. Supplying lower level parameters will permit tweaks to the higher level functions as well.
}
\usage{
phewasManhattan(d, annotate.phenotype.description=T, ...)

phenotypeManhattan(d, suggestive.line = 0.05, significant.line,
  OR.size = F, OR.direction = F, sizes, annotate.level, 
  y.axis.interval = 5, y.axis.label = expression(-log[10](italic(p))), max.y,
  ...)
  
phenotypePlot(d, max.y,max.x, suggestive.line, significant.line, 
           size.x.labels=9, size.y.labels=9, switch.axis=F, 
           sort.by.value=F, sort.by.category.value=F,
           base.labels=F, annotate.phenotype.description,
           annotate.angle=0, annotate.size=5, annotate.level,
           annotate.phenotype=F, annotate.snp.w.phenotype=F,
           annotate.snp=F, annotate.snp.angle=0,
           annotate.list, annotate.only.largest=T,
           lc.labels=F, x.group.labels=T, x.phenotype.labels=F,
           sizes=F, direction=F, point.size=3,
           use.color=T, color.palette,
           title= paste0("Phenotype Plot ", date()),
           x.axis.label="Phenotypes",
           y.axis.label="Values",
           y.axis.interval=5)
}

\arguments{
  \item{d}{
Data frame containing \code{phenotype} and \code{p}, or in the case of \code{phenotypePlot}, \code{phenotype} and \code{value}.
}
  \item{\dots}{
Further parameters to be passed down to \code{phenotypeManhattan} and \code{phenotypePlot}.
}
  \item{suggestive.line}{
For \code{phenotypeManhattan}, Draw a blue line at this p-value to show a suggestive significance threshold. Defaults to 0.05. A missing or NA value as applicable will result in no line.
For \code{phenotypePlot}, Blue line drawn at \code{value} to emphasize interesting data points. Missing or an NA value will result in no line.
}
  \item{significant.line}{
For \code{phenotypeManhattan}, draw a red line at this p-value to show an adjusted significance threshold. Defaults to suggestive.line divided by the number of non-NA p-values. An NA value will result in no line.
For \code{phenotypePlot}, the red line is drawn at the specified \code{value} to emphasize interesting data points. Missing or an NA value will result in no line.
}
  \item{OR.size}{
Adjust point size based on odds ratios? Requires d to contain the \code{OR} column. Default is FALSE. Mutually exclusive with \code{sizes}. 
}
  \item{OR.direction}{
Adjust point shape based on odds ratio direction? Requires d to contain the \code{OR} column. Default is FALSE.
}
  \item{annotate.level}{
For \code{phenotypeManhattan} or \code{phenotypeManhattan}, annotate points at or below this \code{p} with a default of \code{significant.line}. Set to \code{NA} to disable adding these annotations.
For \code{phenotypePlot}, minimum \code{value} to annotate points. Default of no annotation with a missing arugment, or set to \code{NA} to disable.
}
  \item{y.axis.interval}{
The interval for y axis labeling. Defaults to 5. Scale is in -log10 units for the Manhattan plots and is the same as \code{value} for \code{phenotypePlot}.
}
  \item{x.axis.label}{
The label for the x axis. Defaults to "Phenotypes".
}
  \item{y.axis.label}{
The label for the y axis. Defaults to a stylized "-log10(p)" for the Manhattan plots and "Values" for \code{phenotypePlot}.
}
  \item{max.y}{
For the Manhattan plots, the maximum \code{-log10(p)} on the axis. Defaults to the larger of 4.4 or the maximum \code{-log10(p)}  in the supplied data.
For \code{phenotypePlot}, the maximum \code{value} for the plot. If missing, defaults to the maximum value in the supplied data.
}
  \item{max.x}{
Maximum number of \code{phenotype}s for the plot. If missing, defaults to the total number.
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
Sort the plot by highest to lowest \code{-log10(p)} or \code{value}? Defaults to false.
}
  \item{sort.by.category.value}{
Sort the plot by \code{-log10(p)} or \code{value} within each phenotype group? Requires group attributes in \code{d}. Defaults to false.
}
  \item{base.labels}{
Use base ggplot2 labels? Defaults to FALSE, which uses ggrepel.
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
Adjust point size based on \code{size} in \code{d}? Default is FALSE. Mutually exclusive with \code{OR.size} available in the Manhattan plots. 
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
Title for the plot. Defaults to "Phenotype Plot" with the date and time.
}
}
\details{
These functions are nested: \code{phewasManhattan} calls \code{phenotypeManhattan} which calls \code{phenotypePlot}.

\code{phewasManhattan} generates a PheWAS Manhattan plot and can be used with a result data frame from \code{phewas} using phecodes. This function only adds PheWAS code groupings and descriptions as applicable before calling \code{phenotypeManhattan}.

\code{phenotypeManhattan} performs the transformations to the data and paramters for a -log10 scale. It is useful in cases of \code{phewas} results using other phenotypes besides phecodes.

\code{phenotypePlot} does the actual plotting and can be used to visualize non-p-value data. 

All of these functions return \code{\link[ggplot2:ggplot2-package]{ggplot2}} objects. Aspects of the plots can be modified post hoc, for instance to adjust the point size scale one can use \code{plot + scale_size(<parameters>)}
}
\value{
A \code{ggplot2} plot. It contains PheWAS codes against their -log10 transformed p-values for \code{phewasManhattan}, generic phenotypes versus -log10 transformed p-values for \code{phenotypeManhattan}, and generic phenotypes versus generic values for \code{phenotypePlot}.
}

\seealso{
\code{\link[PheWAS:phewas-package]{The PheWAS package}} for a complete example.
}
\keyword{ hplot }

