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
#' phewasManhattan(PheWAS:::test_phewas)
phewasManhattan <-
  function(d, annotate.phenotype.description=T, ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) 
      stop("Data input must contain columns phenotype and p.")
    if(!is(d$phenotype, 'character')) {

      
        stop("Non-character phenotypes passed in, so an accurate
             phecode mapping is not possible.")
      
    }
    #Check to see if it looks 0-padded
    if(min(nchar(d$phenotype))<3) stop('Phenotypes with length <3 observed,ensure they are are 0-padded')

    #Add the groups and phecode descriptions as requested

    d=addPhecodeInfo(d,groupnums =T, groupcolors = T)

    phenotypeManhattan(d, annotate.phenotype.description=annotate.phenotype.description, ...)
  }
