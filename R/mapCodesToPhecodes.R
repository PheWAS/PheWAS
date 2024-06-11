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
#' \code{\link[PheWASmaps:phecode_map]{PheWASmaps::phecode_map}} supports ICD9CM
#' (map v1.2) and ICD10CM (map 2018 beta). If \code{NULL}, it will skip the
#'  mapping codes to phecodes step. This may be useful if one is seeking to
#'   expand or roll up an existing set of phecodes.
#' @param rollup.map Data frame with columns \code{code}, and
#' \code{phecode_unrolled}. Each row represents a mapping from a specific
#'  phecode to all parent phecodes. The default map
#'  \code{\link[PheWASmaps:phecode_rollup_map]{PheWASmaps::phecode_rollup_map}} is the
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
#' @import data.table
#' @examples diabetes_billing=data.frame(id=1:3,vocabulary_id=
#' c("ICD9CM","ICD9CM","ICD10CM"),code=c("250.00","250.01","E11.00"))
#' phecodes=mapCodesToPhecodes(diabetes_billing)
mapCodesToPhecodes <-
  function(input,
           vocabulary.map=PheWASmaps::phecode_map,
           rollup.map=PheWASmaps::phecode_rollup_map,
           make.distinct=TRUE) {
    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[["code"]]) %in% c("character")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}

    if(!is.null(vocabulary.map)){
       if(!is.data.table(input)){
      withCallingHandlers(output <- inner_join(input,vocabulary.map,by=c("vocabulary_id","code")),
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
  #Remove old columns
        }else if(is.data.table(input)){
          setkey(input, vocabulary_id, code)
         vocabulary.map <- as.data.table(vocabulary.map)
          setkey(vocabulary.map, vocabulary_id, code)
        withCallingHandlers(output <- input[vocabulary.map, on = .(vocabulary_id, code), nomatch = NULL, allow.cartesian = TRUE],
                            warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
        }
      output <- output %>% select(-code,-vocabulary_id) %>% rename(code=phecode)
    } else {
      #Warn if the vocabulary IDs are not phecodes
      if(sum(input$vocabulary_id!="phecode")!=0) {stop("Phecode mapping was not requested, but the vocabulary_id of all codes is not 'phecode'")}
      #Prepare for just the phecode expansion
      output=input %>% filter(vocabulary_id=="phecode") %>% select(-vocabulary_id)
    }
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    #Perform the rollup
    if(!is.null(rollup.map)) {
     if(!is.data.table(output)){
      withCallingHandlers(output <- inner_join(output ,rollup.map,by="code"),
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
     } else if(is.data.table(output)){
       setkey(output, 'code')
        #setkey(rollup.map, code)
        withCallingHandlers(output <- output[rollup.map, on = .(code), nomatch = NULL, allow.cartesian = TRUE],
                          warning = function(w) { if (grepl("coercing into character vector", w$message)) {invokeRestart("muffleWarning")}})
    }
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
