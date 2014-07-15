generateExample <- function(n=5000,phenotypes.per=10,hit="335") {
  if(!require(MASS)) stop("Package MASS required to generate artificial association")
  phewas_code=unique(phemap$phewas_code)
  #Exclude the code to add
  phewas_code=phewas_code[phewas_code!=hit]
  #Assign individuals random phenotypes
  random=data.frame(id=rep.int(1:n,phenotypes.per),phewas_code="",count=0)
  random$phewas_code=sample(phewas_code,nrow(random),replace=TRUE)
  #Create the signal
  signal=as.data.frame(mvrnorm(n=n,mu=c(0,0),Sigma=rbind(c(.5,.1),c(.1,1))))
  names(signal)=c("phenotype","rsEXAMPLE")
  #Normalize the phenotype to logical and the genotype to 0,1,2.
  signal$phenotype=signal$phenotype>.25
  signal$rsEXAMPLE=signal$rsEXAMPLE+abs(min(signal$rsEXAMPLE))
  signal$rsEXAMPLE=floor(signal$rsEXAMPLE*2.99/max(signal$rsEXAMPLE))
  signal$id=1:n
  signal$count=0
  signal$phewas_code=hit
  random=rbind(random,signal[signal$phenotype,c("id","phewas_code","count")])
  random=merge(random,phemap)
  random$count= rpois(nrow(random),4)
  random[random$phewas_code==hit,]$count=random[random$phewas_code==hit,]$count+2
  random=random[random$count>0,]
  return(list(id.icd9.count=random[,c("id","icd9","count")],genotypes=signal[,c("id","rsEXAMPLE")]))
}