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
Version: \tab 0.99.5-2\cr
Date: \tab 2019-08-15\cr
License: \tab GPL3\cr
}
This package contains the mappings and information for phecodes version 1.2 (ICD9CM) and the 2018 beta version of ICD10CM.
With the package installed, the command \code{vignette("PheWAS-package")} will display the vignette.
The examples on the PheWASExamples github repository provide some helpful guides, especially \link[https://github.com/PheWAS/PheWASExamples/tree/master/BMI FTO - v1.0]{the V1.0 walkthrough here}.
\itemize{
\item \code{\link[PheWAS:mapCodesToPhecodes]{mapCodesToPhecodes}} allows users to transform their ICD9CM and ICD10CM codes into phecodes.
\item \code{\link[PheWAS:addPhecodeInfo]{addPhecodeInfo}} allows users to add the description for phecodes.
\item \code{\link[PheWAS:createPhenotypes]{createPhenotypes}} takes code count data produces output suitable for PheWAS anaylsis.
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
#Extract the three parts from the returned list
id.vocab.code.count=ex$id.vocab.code.count
genotypes=ex$genotypes
id.sex=ex$id.sex
#Create the phecode table- translates the codes, adds 
#exclusions, and reshapes to a wide format.
#Sum up the counts in the data where applicable.
phenotypes=createPhenotypes(id.vocab.code.count, 
  aggregate.fun=sum, id.sex=id.sex)
#Combine the data
data=inner_join(inner_join(phenotypes,genotypes),id.sex)
#Run the PheWAS
results=phewas_ext(data,phenotypes=names(phenotypes)[-1],genotypes=c("rsEXAMPLE"),
  covariates=c("sex"), cores=1)
#Plot the results
phewasManhattan(results, 
  title="My Example PheWAS Manhattan Plot")
#Add PheWAS descriptions
results_d=addPhecodeInfo(results)
#List the top 10 results
results_d[order(results_d$p)[1:10],]
#Create a nice interactive table (eg, in RStudio)
phewasDT(results)
}
}

