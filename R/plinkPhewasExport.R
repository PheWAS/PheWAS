plinkPhewasExport <- function(phenotypes, file="plink.pheno", translateIDs=TRUE) {
  if(translateIDs) {
    #Create the new ID columns and remove the old one (whatever it was called)
    phenotypes$FID=phenotypes[[1]]
    phenotypes$IID=phenotypes[[1]]
    phenotypes=phenotypes[,-1]
  }
  #Make sure the IDs are in the proper order
  pheno.out=phenotypes %>% select(FID, IID, everything())
  pheno.out=lapply(pheno.out,function(x){
    if(class(x)=="logical"){
      x=ifelse(is.na(x),-9,ifelse(x,2,1))
    }
    x
  })
  write.table(data.frame(pheno.out,check.names = F), file = file, row.names = FALSE, quote = FALSE, sep="\t")
}