mapICD9toPheWAS <-
function(..., icd9s, add.info=F, keep.icd9=F) {
  if(!missing(icd9s)) {
    if(missing(...)){
      input=tbl_df(data.frame(icd9=icd9s))
    } else {
      input=tbl_df(data.frame(list(...),icd9=icd9s))
    }
  }
  else {
    input=tbl_df(data.frame(...))
  }
  if(sum(names(input)=="icd9")==0) {
    stop("Must supply a data frame with an 'icd9' column or specify the icd9s parameter")
  }
  if(class(input$icd9) %in% c("integer","numeric")) {stop("Numeric ICD-9 codes passed in, so an accurate mapping is not possible. E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
  #merge the tables
  output = inner_join(input,phemap,by="icd9")
  if(add.info){
    output=inner_join(output,pheinfo,by="phewas_code")
  }
  if(!keep.icd9) {
    output=output %>% select(-icd9)
  }
  output
}
