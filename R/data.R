
#' Sample Data
#'
#'Sample data for documentation and testing
#'
#' @format ## `sample_data`
#' list of 3
#' \describe{
#' \item{id.vocab.code.count}{A data frame containing personID, vocabulary_ID, code, and count}
#' \item{genotypes}{a dataframe containing the person_ID and rsEXAMPLE status}
#' \item{id.sex}{a data frame containing the person_ID and their gender}
#' ...
#' }
"sample_data"

#' example_phewas
#' 
#' Example Data used in documentation
#'
#' @format ## `example_phewas`
#' data frame with 1825 rows and 17 columns
#' \describe{
#' \item{phenotype}{Sample Phenotypes}
#' \item{snp}{sample SNP}
#' \item{covariates}{Sample Covariates}
#' \item{beta}{the beta of the sample analysis}
#' \item{SE}{Standard Error of the sample analysis}
#' \item{OR}{Odds Ratio of the sample analysis}
#' \item{p}{P value of the sample analysis}
#' \item{type}{Logistic Regression}
#' \item{n_total}{number of samples in sample analysis}
#' \item{n_cases}{number of cases in sample analysis}
#' \item{n_controls}{number of controls in sample analysis}
#' \item{HWE_p}{Hardy-Weinberg}
#' \item{allele_freq}{Allele Frequency}
#' \item{n_no_snp}{Fill In Later}
#' \item{formula}{Fill In Later}
#' \item{expanded_formula_note}{Fill in Later}
#' ...
#' }
"example_phewas"