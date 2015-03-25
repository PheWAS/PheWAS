\name{phewasMetaModels}
\alias{phewasMetaModels}
\title{
Perform meta-analysis of PheWAS results
}
\description{
This function wraps the meta package's \code{\link[meta]{metagen}} function to meta analyze \code{phewas} results.
}
\usage{
phewasMetaModels(results, ...)
}
\arguments{
  \item{results}{
Data frame containing \code{phewas} results. This data frame should include an additional column, \code{study}, identifying which study the results are for. See details and examples for more information.
}
  \item{\dots}{
Additional parameters to be passed to \code{metagen}.
}
}
\details{
This function uses \code{by} to split \code{results} into groups of studies based on \code{phenotype, snp, } and \code{adjustment}. \code{phewasMeta} forces \code{NA} adjustment values to be a character string "NA" due to \code{by} restrictions on the \code{INDICES}.
Note that \code{update.meta} can be used to alter individual \code{metagen} objects.
}
\value{
A \code{by} object with named entries containing the meta analysis models. Each entry is named in the form `phenotype snp adjustment`.
}
\seealso{
\code{\link[PheWAS:phewasMeta]{phewasMeta}} will return a data frame of results similar to \code{phewas} output.
\code{\link[meta:update.meta]{update.meta}} will allow users to alter returned \code{metagen} objects.
}

\examples{
\donttest{
#Generate some example data
ex=generateExample(hit="335")
#Extract the two parts from the returned list
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
#Create the PheWAS code table- translates the icd9s, adds exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count)
#Run the PheWAS
results.1=phewas(phenotypes,genotypes,cores=2,significance.threshold=c("bonferroni"))
#Set up a study identifier
results.1$study="335"
#Perform another PheWAS
ex=generateExample(hit="250.2")
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
phenotypes=createPhewasTable(id.icd9.count)
results.2=phewas(phenotypes,genotypes,cores=2,significance.threshold=c("bonferroni"))
results.2$study="250.2"
#Join the two sets of PheWAS results
results=rbind(results.1,results.2)
#Perform the meta analysis, and do not assume fixed effects.
results.meta=phewasMetaModels(results)
results.meta$`250.2 rsEXAMPLE NA`
}
}
\keyword{ htest }

