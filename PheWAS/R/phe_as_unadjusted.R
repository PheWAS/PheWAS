phe_as_unadjusted <-
  function(phe_gen) {
    #Retrieve the targets for this loop
    phe_o=phe_gen[1]
    phe=phe_o
    gen=phe_gen[2]
    d=data[,c(gen,phe)]
    d=na.omit(d)
    
    assign("last.warning", NULL, envir = baseenv())

    if(length(unique(d[,phe]))<=10 & length(unique(d[,gen]))<=10) {
      type = "chi-square"
      p=chisq.test(table(d[,gen],d[,phe]))$p.value
    } else if(length(unique(d[,phe]))==2) {
      p=t.test(as.formula(paste0(phe," ~ ", gen, collapse="")),data=d)$p.value
      type = "t-test"
    } else {
      stop("Unadjusted models were forced, but the outcome measure was not categorical")
    }
    output=data.frame(phenotype=phe_o,snp=gen,
                      p=p, type=type,
                      note=ifelse(!is.null(warnings()),warnings(),""), stringsAsFactors=F)
    
    
    #Return this to the loop to be merged.
    output
  }
