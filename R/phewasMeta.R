phewasMeta <- function(results, fixed=T, keep.both=T, cores=1, ...) {
  #If parallel, run the parallel version.
  if(cores>1) {
    #Prep the cluster
    cluster=multidplyr::new_cluster(cores)
    #Iterate across all phenotype, snp, adjustment combinations.
    out = results %>% group_by(phenotype,snp,adjustment) %>%  multidplyr::partition(cluster) %>% do(phewas_meta_logic(., fixed=fixed, ...)) %>% collect()
  } else {
    #Otherwise, just use do.
    #Iterate across all phenotype, snp, adjustment combinations.
    out = results %>% group_by(phenotype,snp,adjustment) %>% do(phewas_meta_logic(., fixed=fixed, ...))    
  }
  out
}