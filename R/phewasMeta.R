phewasMeta <- function(results, fixed=T, keep.both=T, cores=1, ...) {
  #Replace NA adjustment values with "NA"- necessary for by
  #results$adjustment = ifelse(is.na(results$adjustment),"NA",results$adjustment)
  
  #If parallel, run the parallel version.
  if(cores>1) {
    if(!require(multidplyr)) {stop("Package 'multidplyr' required for parallelization. Please see https://github.com/hadley/multidplyr for details.")}
    #Prep the cluster
    cluster=parallel::makePSOCKcluster(cores)
    #Iterate across all phenotype, snp, adjustment combinations.
    out = results %>% partition(phenotype,snp,adjustment,cluster=cluster) %>% do(PheWAS:::phewas_meta_logic(., fixed=fixed, ...)) %>% collect()
  } else {
    #Otherwise, just use do.
    #Iterate across all phenotype, snp, adjustment combinations.
    out = results %>% group_by(phenotype,snp,adjustment) %>% do(phewas_meta_logic(., fixed=fixed, ...))
    
  }
}