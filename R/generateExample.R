generateExample <- function(n=5000,phenotypes.per=10,hit="335") {
    phecode=unique(pheinfo$phecode)
    #Exclude the code to add
    phecode=phecode[phecode!=hit]
    #Assign individuals random phenotypes
    random=data.frame(id=rep.int(1:n,phenotypes.per),phecode="",count=0)
    random$phecode=sample(phecode,nrow(random),replace=TRUE)
    #Create the signal
    signal=as.data.frame(MASS::mvrnorm(n=n,mu=c(0,0),Sigma=rbind(c(.5,.1),c(.1,1))))
    names(signal)=c("phenotype","rsEXAMPLE")
    #Normalize the phenotype to logical and the genotype to 0,1,2.
    signal$phenotype=signal$phenotype>.25
    signal$rsEXAMPLE=signal$rsEXAMPLE+abs(min(signal$rsEXAMPLE))
    signal$rsEXAMPLE=floor(signal$rsEXAMPLE*2.99/max(signal$rsEXAMPLE))
    signal$id=1:n
    signal$count=0
    signal$phecode=hit
    random=rbind(random,signal[signal$phenotype,c("id","phecode","count")])
    random=merge(random,phemap)
    random$count= rpois(nrow(random),4)
    random[random$phecode==hit,]$count=random[random$phecode==hit,]$count+2
    random=random[random$count>0,]
    return(list(id.icd9.count=random[,c("id","icd9","count")],genotypes=signal[,c("id","rsEXAMPLE")]))
}