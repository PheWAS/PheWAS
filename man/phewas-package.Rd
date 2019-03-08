\name{PheWAS-package}
\alias{PheWAS-package}
\alias{PheWAS}
\docType{package}
\title{
PheWAS analysis methods
}
\description{
This package provides the tools necessary to perform a PheWAS analysis.
}
\details{
\tabular{ll}{
Package: \tab PheWAS\cr
Type: \tab Package\cr
Version: \tab 0.99.3-1\cr
Date: \tab 2019-03-08\cr
License: \tab GPL3\cr
}
This package contains the mappings and information for phecodes version 1.2 (ICD9CM) and the 2018 beta version of ICD10CM.
With the package installed, the command \code{vignette("PheWAS-package")} will display the vignette
\itemize{
\item \code{\link[PheWAS:mapCodesToPhecodes]{mapCodesToPhecodes}} allows users to transform their ICD9CM and ICD10CM codes into phecodes.
\item \code{\link[PheWAS:addPhecodeInfo]{addPhecodeInfo}} allows users to add the description for phecodes.
\item \code{\link[PheWAS:createPhewasTable]{createPhewasTable}} takes ICD9 count data produces output suitable for PheWAS anaylsis.
\item \code{\link[PheWAS:phewas]{phewas}} allows users to perform a PheWAS analysis with genetics, ICD9 codes, phecodes, etc.
\item \code{\link[PheWAS:PheWAS_Plotting]{phewasManhattan}} is a high level plot function that will generate a ggplot2 of PheWAS results.
\item \code{\link[PheWAS:PheWAS_Plotting]{phenotypeManhattan}} provides a more general interface for plotting phenotype based Manhattan plots.
\item \code{\link[PheWAS:PheWAS_Plotting]{phenotypePlot}} provides a complex interface for plotting many types of phenotype data.
\item \code{\link[PheWAS:phewasMeta]{phewasMeta}} will perform a meta analysis of \code{phewas} results.
\item \code{\link[PheWAS:generateExample]{generateExample}} creates data for example or testing purposes.
}
}
\author{
Robert Carroll <robert.carroll@vumc.org>
}
\references{
The original PheWAS manuscript can be found here: \url{http://www.ncbi.nlm.nih.gov/pubmed/20335276}
}
\keyword{ package}
\keyword{ utilities }
\seealso{
\code{\link[PheWAS:mapICD9ToPhecodes]{mapICD9ToPhecodes}}

\code{\link[PheWAS:addPhecodeInfo]{addPhecodeInfo}}

\code{\link[PheWAS:createPhewasTable]{createPhewasTable}}

\code{\link[PheWAS:phewas]{phewas}}

\link[PheWAS:PheWAS_Plotting]{PheWAS plotting methods}

\code{\link[PheWAS:phewasMeta]{phewasMeta}}
}
\examples{
\donttest{
#Load the PheWAS package
library(PheWAS)
#Set the random seed so it is replicable
set.seed(1)
#Generate some example data
ex=generateExample()
#Extract the two parts from the returned list
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
#Create the phecode table- translates the icd9s, adds 
#exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count)
#Run the PheWAS
results=phewas(phenotypes,genotypes,cores=1,
  significance.threshold=c("bonferroni"))
#Plot the results
phewasManhattan(results, annotate.angle=0,
  title="My Example PheWAS Manhattan Plot")
#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]
}
}

