#' Map codes to phecodes
#'
#' This function takes a data frame with codes and maps them to phecodes. It can
#'  support aribitrary maps of codes to phecodes and phecode rollup maps. The
#'  included mappings are from ICD9CM and ICD10CM.
#'
#' @param input Data frame containing \code{vocabulary_id} and \code{code}
#' columns. These columns specify the vocabulary used in each row, eg ICD9CM or
#'  ICD10CM, and the code to be translated. \code{code} must be a character or
#'  factor to ensure proper conversion (ICD9CM codes lose specificity as numeric).
#' @param vocabulary.map Data frame with columns \code{vocabulary_id},
#' \code{code}, and \code{phecode}. Each row represents a mapping from a
#' specific vocabulary and code to a specific phecode. The default map
#' \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}} supports ICD9CM
#' (map v1.2) and ICD10CM (map 2018 beta). If \code{NULL}, it will skip the
#'  mapping codes to phecodes step. This may be useful if one is seeking to
#'   expand or roll up an existing set of phecodes.
#' @param rollup.map Data frame with columns \code{code}, and
#' \code{phecode_unrolled}. Each row represents a mapping from a specific
#'  phecode to all parent phecodes. The default map
#'  \code{\link[PheWAS:phecode_rollup_map]{PheWAS::phecode_rollup_map}} is the
#'   complete rollup map for phecode map v1.2. If \code{NULL}, it will skip the
#'   rollup step. This may be useful if one is seeking to only consider the
#'   directly mapped phecodes.
#' @param make.distinct Boolean value. Should duplicate rows be removed during
#' mapping? Default is \code{TRUE}. Useful to reduce data size, especially when
#'  another column, eg date of code, is provided.
#'
#' @return A data frame containing the columns in \code{input}, except the
#' original \code{code} and \code{vocabulary_id} columns have been replaced
#' with \code{phecode} now containing the phecode as character.
#' @export
#'
#' @examples diabetes_billing=data.frame(id=1:3,vocabulary_id=
#' c("ICD9CM","ICD9CM","ICD10CM"),code=c("250.00","250.01","E11.00"))
#' phecodes=mapCodesToPhecodes(diabetes_billing)
#' phecodes
#' addPhecodeInfo(phecodes)
#'
mapCodesToPhecodes <-
  function(input,
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           make.distinct=TRUE) {
    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[["code"]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}

    if(!is.null(vocabulary.map)){
      #Perform the direct map
      withCallingHandlers(output <- inner_join(input,vocabulary.map,by=c("vocabulary_id","code")),
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
      #Remove old columns
      output = output %>% select(-code,-vocabulary_id) %>% rename(code=phecode)
    } else {
      #Warn if the vocabulary IDs are not phecodes
      if(sum(input$vocabulary_id!="phecode")!=0) {warning("Phecode mapping was not requested, but the vocabulary_id of all codes is not 'phecode'")}
      #Prepare for just the phecode expansion
      output=input %>% filter(vocabulary_id=="phecode") %>% select(-vocabulary_id)
    }
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    #Perform the rollup
    if(!is.null(rollup.map)) {
      withCallingHandlers(output <- inner_join(output ,rollup.map,by="code"),
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
      output = output %>% select(-code) %>% rename(phecode=phecode_unrolled)
      #Make distinct
      if(make.distinct) {output = distinct(output)}
    } else {
      #Rename output column to phecode
      output = output %>% rename(phecode=code)
    }

    #Return the output
    output
  }
