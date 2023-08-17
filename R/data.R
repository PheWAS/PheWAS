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

