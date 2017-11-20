phewasMeta <- function(results, fixed=T, keep.both=T, ...) {
  #Replace NA adjustment values with "NA"- necessary for by
  results$adjustment=ifelse(is.na(results$adjustment),"NA",results$adjustment)
  #Iterate across all phenotype, snp, adjustment combinations.
  out= results %>% group_by(phenotype,snp,adjustment) %>% do(phewas_meta_logic(., ...))
  #bind_rows(out)
}