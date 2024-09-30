#' Get Pheno is a helper function for the long table PheWAS format
#'
#' @param pheno Phecode Counts table
#' @param demos Demos table - mainly requires sex
#' @param phecode Phecode info table
#' @param MCC Minimum code count, default is 2
#'
#' @return a data table containing the demographic data with an additional column
#' named 'pheno' where 1 is case, 0 is control, -9 is excluded


get_pheno <- function(pheno, demos, phecode, MCC=2){
  
  cur_phecode <- phecode
  
  p <- pheno[phecode == cur_phecode, c('GRID','N')]
  p <- data.table(p, key = "GRID")
  # default pheno to -9 (flags N<MCC)
  p[, pheno := -9]
  p[N >= MCC, pheno := 1]
  p[, N := NULL]#delete column N
  
  p=merge(demos, p,by="GRID",all.x=TRUE)
  p[is.na(pheno), pheno := 0] #set pheno to 0 for controls
  p[pheno == -9, pheno := NA] #set pheno to NA (not a case or control)
  
  return(p)
}