phewasDT <- function(results, addPhecodeInfo=T, OR.digits=2, p.digits=3, ...){
  if(addPhecodeInfo==T) results = results %>% addPhecodeInfo()
  dt.res=results %>% filter(!is.na(p)) %>% 
    transmute(phenotype,description,snp,OR=round(OR,OR.digits),p=signif(p,p.digits),n_total,n_cases,n_controls)
  DT::datatable(dt.res, rownames=F,...)
}