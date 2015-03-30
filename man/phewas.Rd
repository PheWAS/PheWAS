\name{phewas}
\alias{phewas}

\title{
Function to perform a PheWAS analysis
}
\description{
This function will perform a PheWAS analysis, optionally adjusting for other variables. It is parallelized using snowfall.
}
\usage{
phewas(phenotypes, genotypes, data, covariates = c(NA), adjustments = list(NA), 
  outcomes, predictors, cores = 1, additive.genotypes = T, 
  significance.threshold, alpha = 0.05, unadjusted = F, return.models = F,
  min.records = 20, MASS.confint.level=NA, quick.confint.level)
}
\arguments{
  \item{phenotypes}{
The names of the outcome variables in data under study. These can be logical (for logistic regression) or continuous (for linear regression) columns. Can alternatively be a data frame of phenotypes, see Details for more information.
}
  \item{genotypes}{
The names of the prediction variables in data under study. These can be logical or continuous. Can alternatively be a data frame of genotypes, see Details for more information.
}
  \item{data}{
Data frame containing all variables for the anaylsis. Omitted if using data frames for other parameters, see Details for more information
}
  \item{covariates}{
The names of the covariates to appear in every analysis. Can alternatively be a data frame of covariates, see Details for more information.
}
  \item{adjustments}{
A list containing special one-off adjustments for the analysis. \code{list(NA)} will yield no special adjustments. Use the \code{covariates} parameter to adjust the analyses under normal circumstances. See Details for more information. 
}
  \item{outcomes}{
An alternate to \code{phenotypes}. It will be ignored if \code{phenotypes} exists.
}
  \item{predictors}{
An alternate to \code{genotypes}. It will be ignored if \code{genotypes} exists.
}
  \item{cores}{
The number of cores to use in the snowfall implementation. If \code{cores=1},  \code{lapply} will be used instead.
}
  \item{additive.genotypes}{
Are additive genotypes being supplied? If so, it will attempt to calculate allele frequencies and HWE values. Default is TRUE.
}
  \item{significance.threshold}{
A vector of desired significance thresholds to calculate. Can include "p-value","bonferroni","fdr","simplem-genotype","simplem-phenotype","simplem-product". Note that simpleM based methods can be time intensive. See details for more information.
}
  \item{alpha}{
The base alpha for significance calculations.
}
  \item{unadjusted}{
Use Chi-Square and t-tests instead? This is a much simpler implementation. Defaults FALSE.
}
  \item{return.models}{
Return a list the complete models, with the names equal to the string formula used to create them, as well as the results. Default is FALSE.
}
  \item{min.records}{
The minimum number of records to perform a test. For logistic regression, there must be at least this number of each cases and controls, for linear regression this total number of records. Default is 20.
}
  \item{MASS.confint.level}{
Uses the \code{MASS} package and the \code{confint} function to calculate a confidence interval at the specified level. \code{MASS::confint} uses a profile likelihood method, which takes some time to compute. Output is stored in the \code{lower} and \code{upper} columns. Logistic models will report OR CIs and linear models will report beta CIs. Default is NA, which does not calculate confidence intervals.
}
  \item{quick.confint.level}{
Calculate a confidence interval based on \code{beta + or - qnorm * SE}. Output is stored in the \code{lower.q} and \code{upper.q} columns. Logistic models will return have the exponentiated OR confidence intervals.
}
}
\details{
The complete data frame can be passed in using the \code{data} parameter with name vectors in \code{phenotypes, genotypes, covariates, and adjustments} parameters. Alternatively, \code{phenotypes, genotypes, covariates, and adjustments} can each be data frames. They will be merged using the shared columns between \code{phenotypes} and \code{genotypes}, ideally being an ID column.

\code{covariates} are those variables that are included in every model, e.g., age and gender. \code{adjustments} are additional variables that can be used to compare models by adjusting for potentially confounding factors. Including \code{NA} in the list will perform a set of analyses with no special adjustments, single names can be used for single adjustment, and a vector of names can be used for multiple special adjustments at once. An \code{adjustments} parameter of \code{list(NA, "BMI", c("BMI", "smoking"))} would adjust for no extra variable in all models, followed by adjusting additionally for BMI, and then perform another round of anaylses using the BMI and smoking status variables. 

These results can be directly plotted using the \code{\link[PheWAS:phewasManhattan]{phewasManhattan}} function, assuming that models are not returned. If they are, the \code{results} item of the returned list needs to be used.
}
\value{
The following are the default rows included in the returned data frame. The attributes of the returned data frame contain additional information about the anaylsis.
\item{phenotype}{The outcome under study}
\item{snp}{The predictor under study}
\item{adjustment}{The one off adjustment used}
\item{beta}{The beta coefficient for the predictor}
\item{SE}{The standard error for the beta coefficient}
\item{lower.p}{The lower bound of the quick confidence interval, if requested}
\item{upper.p}{The upper bound of the quick confidence interval, if requested}
\item{lower}{The lower bound of the \code{MASS::confint} confidence interval, if requested}
\item{upper}{The upper bound of the \code{MASS::confint} confidence interval, if requested}
\item{OR}{For logistic regression, the odds ratio for the predictor}
\item{p}{The p-value for the predictor}
\item{type}{The type of regression model used}
\item{n_total}{The total number of records in the analysis}
\item{n_cases}{The number of cases in the analysis (logical outcome only)}
\item{n_controls}{The number of controls in the analysis (logical outcome only)}
\item{HWE_p}{The Hardy-Weinberg equilibrium p-value for the predictor, assuming 0,1,2 allele coding}
\item{allele_freq}{The allele frequency in the predictor for the coded allele}
\item{n_no_snp}{The number of records with a missing predictor}
\item{note}{Additional warning or error information}
If there are any requested significance thresholds, boolean variables will be included reporting significance.
If \code{return.models=T}, a list is returned. The named item \code{results} contains the above data frame. The named item \code{models} contains a list of the models generated in the analysis. To distinguish models, the list is named by the full formula used in generation.
}

\author{
Robert Carroll
}

\seealso{
\code{\link[PheWAS:createPhewasTable]{createPhewasTable}} 
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
results=phewas(phenotypes,genotypes,cores=4)
}
}
\keyword{ models }

