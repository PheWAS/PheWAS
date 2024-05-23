#' A function to join all the covariate and phenotype files.
#' 
#' This function takes up to three tables.  These tables include thee phenotype
#'  table (the output of createPhenotype), id.sex (optional), and covar 
#'  (data$genotypes in sample data, 
#' but this can be whatever dependent variable and 
#' covariates you're looking for)
#'
#' @param pheno a wide table with columns representing phecodes and rows 
#' representing people
#' @param id.sex an nx2 table with person_ids and gender information
#' @param covar a table of variable width consisting of person_IDs and whatever 
#' other information the user requires as covariates or dependent variabes.
#'
#' @return a wide data table consisting of person_ids as rows and covariates and
#'  phecodes as columns
#' @export
#'
#' @examples
#' phenotype_data <- createPhenotypes(sample_data$id.vocab.code.count, id.sex = 
#' sample_data$id.sex)
#' joinCovar(phenotype_data, sample_data$id.sex, sample_data$genotypes)
joinCovar <- function(pheno, id.sex, covar){
if(missing(pheno)){
  stop('Phenotype table not found')
}
  if(!missing(id.sex)){
  if(!missing(covar)){
    print('ID.Sex and covar present')
  final_data <- dplyr::inner_join(dplyr::inner_join(id.sex, covar),  
                                  pheno)
  }else{
    print('ID.Sex Present')
    final_data <- dplyr::inner_join(pheno, id.sex)
  }
}else if(!missing(covar)){
  print('Covar presenet')
  final_data <- dplyr::inner_join(pheno, covar)
} else {
  print('No inputs to join, returning phenotype file')
  final_data <- pheno
}
  return(final_data)
}