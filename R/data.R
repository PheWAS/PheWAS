#' Gender Restriction Table
#'
#' A data map that connects phecodes to their gender tendencies.
#'
#' @format ## `gender_restriction`
#' data frame with 1817 rows and 3 columns
#' \describe{
#' \item{phecode}{The phecodes}
#' \item{male_only}{a boolean indicating whether the phecode is male only}
#' \item{female_only}{a boolean indicating whether the phecode is female only}
#' ...
#' }
"gender_restriction"

#' Phecode Exclude Table
#'
#' A data map that connects phecodes to their eexclusions
#'
#' @format ## `phecode_exclude`
#' data frame with 30028 rows of 2 columns
#' \describe{
#' \item{code}{the reference code}
#' \item{exclusion_criteria}{the phecode required to exclude the reference code}
#' ...
#' }
"phecode_exclude"

#' Map from ICD codes to Phecodes
#'
#' A data map that connects icd codes to phecodes
#'
#' @format ## `phecode_map`
#' a data frame with 84532 rows of 3 columns
#' \describe{
#' \item{vocabulary_id}{the ICD source of the ICD code, eg. ICD9CM, ICD10CM}
#' \item{code}{The ICD code}
#' \item{phecode}{The applicable phecode}
#' ...
#' }
"phecode_map"

#' Phecode Rollup Map
#'
#'a data map that contains phecode parent-child pairings
#'
#'@format ## `phecode_rollup_map`
#'a data frame with 3319 rows of 2 columns
#'\describe{
#'\item{code}{the reference phecode}
#'\item{phecode_unrolled}{the parent phecode}
#'...
#'}
"phecode_rollup_map"

#' Pheinfo map
#'
#' a data map that contains information on each phecode
#'
#' @format ## `pheinfo`
#' a data frame with 1817 rows of 5 columns
#' \describe{
#' \item{phecode}{the reference phecode}
#' \item{description}{a string describing the phecode}
#' \item{groupnum}{TBD}
#' \item{group}{what type of disease the phecode is}
#' \item{color}{TBD}
#' ...
#' }
"pheinfo"

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


#' phecode_map_icd10
#'
#' Dataframe containing phecode mappings for ICD10
#' 
#' This data frame maps each ICD10 code to its directly mapped phecode(s).
#' Can be provided to \code{\link[PheWAS:createPhenotypes]{createPhenotypes}} 
#' in the \code{vocabulary.map} parameter to convert ICD10 codes for use in 
#' the PheWAS methods.
#' It is the map 1.2b1.
#' \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}} contains mapping for 
#' ICD9CM and ICD10CM (ie, US use cases).
#'
#' @format ## `phecode_map_icd10`
#'  A data frame with 9505 observations on the following 3 variables.
#' \describe{
#' \item{\code{vocabulary_id}}{Character vector representing the vocabulary of 
#' the code- always 'ICD10'}
#' \item{\code{code}}{Character vector representing the specific code to be 
#' mapped}
#' \item{\code{phecode}}{Character vector representing the direct phecode 
#' mapping}
#' ...
#' }
"phecode_map_icd10"
