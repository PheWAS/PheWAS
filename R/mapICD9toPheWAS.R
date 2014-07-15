mapICD9toPheWAS <-
function(..., icd9s, add.info=F, keep.icd9=F) {
  if(!missing(icd9s)) {
    input=data.frame(list(...),icd9=icd9s)
  }
  else {
    input=data.frame(...)
  }
  if(is.null(input$icd9)) {
    stop("Must supply a data frame with an 'icd9' column or specify the icd9s parameter")
  }
  if(class(input$icd9) %in% c("integer","numeric")) {stop("Numeric ICD-9 codes passed in, so an accurate mapping is not possible. E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
  output = merge(input,phemap)
  if(add.info){
    output=merge(output,pheinfo)
  }
  if(!keep.icd9) {
    output=subset(output, select = -icd9)
  }
  output
}
