\name{generateExample}
\alias{generateExample}
\title{
Generate example data for PheWAS
}
\description{
\code{generateExample} creates simulated ICD9CM/ICD10CM, id.sex, and genotype data to test PheWAS methods. Requires the package \code{MASS}.
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
The phecode for which an artificial hit will be generated.
}
}
\details{
The output from this function will be random results with a significant association found between the genotype data and the selected phecode. It does not represent actual treatment data, though some highly connected phenotypes can also have strong associations due to the data creation methods.
}
\value{
List with three named items:
\item{id.vocab.code.count}{A data frame of columns id, vocabulary_id, ICD code, and count suitable for transformation by \code{createPhenotypes}.}
\item{genotypes}{A data frame of columns id and rsEXAMPLE, suitable for use in \code{phewas} after joining with phenotypic data.}
\item{id.sex}{A data frame of columns id and sex, suitable for use in \code{phewas} and \code{createPhenotypes}.}
}

\keyword{ datagen }
