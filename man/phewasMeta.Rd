\name{phewasMeta}
\alias{phewasMeta}
\title{
Perform meta-analysis of PheWAS results
}
\description{
This function wraps the meta package's \code{\link[meta]{metagen}} function to meta analyze \code{phewas} results.
}
\usage{
phewasMeta(results, fixed=T, keep.both=T, ...)
}
\arguments{
  \item{results}{
Data frame containing \code{phewas} results. This data frame should include an additional column, \code{study}, identifying which study the results are for. See details and examples for more information.
}
  \item{fixed}{
Use the fixed effects model results in columns \code{p, OR, beta,} and \code{SE}? Default is \code{TRUE}.
}
  \item{keep.both}{
If \code{FALSE}, remove the specified fixed effects and random effects results columns (e.g., \code{p.fixed} and \code{p.random}). Default is \code{TRUE}.
}
  \item{\dots}{
Additional parameters to be passed to \code{metagen}.
}
}
\details{
This function uses \code{by} to split \code{results} into groups of studies based on \code{phenotype, snp, } and \code{adjustment}. \code{phewasMeta} forces \code{NA} adjustment values to be a character string "NA" due to \code{by} restrictions on the \code{INDICES}.
The results dataframe must contain an additional column, \code{study}, which contains a study identifier for the meta analysis. Adding a \code{study} column for each set of phewas results followed by an \code{rbind} joining them together would create a functional \code{results} parameter. See the examples section.
}
\value{
A data frame containing meta-analysis results. Includes most phewas result columns and summary information on the meta analysis.
}
\seealso{
\code{\link[PheWAS:phewasMetaModels]{phewasMetaModels}} will return a \code{by} object of complete \code{metagen} objects.
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
results.1=phewas(phenotypes,genotypes,cores=4,significance.threshold=c("bonferroni"))
#Set up a study identifier
results.1$study="335"
#Perform another PheWAS
ex=generateExample(hit="250.2")
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
phenotypes=createPhewasTable(id.icd9.count)
results.2=phewas(phenotypes,genotypes,cores=4,significance.threshold=c("bonferroni"))
results.2$study="250.2"
#Join the two sets of PheWAS results
results=rbind(results.1,results.2)
#Perform the meta analysis, and do not assume fixed effects.
results.meta=phewasMeta(results, fixed=FALSE, keep.both=FALSE)
}
}
\keyword{ htest }

