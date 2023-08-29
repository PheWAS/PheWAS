#' Fill in later
#'
#' @param d n/a
#' @param annotate.phenotype.description n/a
#' @param ... n/a
#'
#' @return n/a
#' @export
#'
#' @examples 
#' 
#' data <- sample_data
#'phenotype_data <- createPhenotypes(data$id.vocab.code.count, 
#'id.sex = data$id.sex) 
#'final_data <- dplyr::inner_join(dplyr::inner_join(
#'data$id.sex, data$genotypes),  phenotype_data)
#'test_phewas <- phewas_ext(names(phenotype_data)[-1], 
#'                          genotypes = c('rsEXAMPLE'), covariates = 'sex', 
#'                          data = final_data)
#' phewasManhattan(test_phewas)
phewasManhattan <-
  function(d, annotate.phenotype.description=T, ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) 
      stop("Data input must contain columns phenotype and p.")
    if(!is(d$phenotype, 'character')) {

      if(is(d$phenotype, 'factor')) {
        warning("Factor phenotype input mapped to characters")
        d$phenotype=as.character(d$phenotype)
      } else {
        stop("Non-character or non-factor phenotypes passed in, so an accurate
             phecode mapping is not possible.")
      }
    }
    #Check to see if it looks 0-padded
    if(min(nchar(d$phenotype))<3) warning("Phenotypes with length <3 observed,
                                          ensure they are are 0-padded
                                          (e.g., \"008\")")

    #Add the groups and phecode descriptions as requested

    d=addPhecodeInfo(d,groupnums =T, groupcolors = T)

    phenotypeManhattan(d, annotate.phenotype.description=annotate.phenotype.description, ...)
  }
