contr.phewas=function(x){
  y=contrasts(x)
  colnames(y)=paste0('-',colnames(y))
  y
}