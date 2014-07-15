\name{generateExample}
\alias{generateExample}
\title{
Generate example data for PheWAS
}
\description{
\code{generateExample} creates simulated ICD9 and genotype data to test PheWAS methods. Requires the package \code{MASS}.
}
\usage{
generateExample(n=5000,phenotypes.per=10,hit="335")
}
\arguments{
  \item{n}{
The number of individuals in the study.
}
  \item{phenotypes.per}{
The base number of phenotypes per individual.
}
  \item{hit}{
The PheWAS code for which an artificial hit will be generated.
}
}
\details{
The output from this function will be random results with a significant association found between the genotype data and the selected PheWAS code. It does not represent actual treatment data, though some highly connected phenotypes can also have strong associations due to the data creation methods.
}
\value{
List with two named items:
\item{id.icd9.count}{A data frame of columns id, icd9, and count suitable for transformation by \code{createPhewasTable}.}
\item{genotypes}{A data frame of columns id and rsEXAMPLE, suitable for use in \code{phewas}}
}

\keyword{ datagen }
