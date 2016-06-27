plinkPhewasExport <- function(phenotypes, file="plink.pheno", translateIDs=TRUE) {
  if(translateIDs) {
    name_id=names(phenotypes)[1]
    #Create the new ID columns
    phenotypes$FID=phenotypes[[1]]
    phenotypes$IID=phenotypes[[1]]
    #Remove the old id column as long as it wasn't one we need
    if(!(name_id %in% c("FID","IID"))) phenotypes=phenotypes[,-1]
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