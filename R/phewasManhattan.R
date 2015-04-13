phewasManhattan <-
  function(d, add.phewas.descriptions=T, ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) stop("Data input must contain columns phenotype and p.")
    if(class(d$phenotype)!="character") {
      if(class(d$phenotype)=="factor") {
          warning("Factor phenotype input mapped to characters")
          d$phenotype=as.character(d$phenotype)
        } else {
          stop("Non-character or non-factor phenotypes passed in, so an accurate phewas code mapping is not possible.")
        }
    }
    #Check to see if it looks 0-padded
    if(min(nchar(d$phenotype))<3) warning("Phenotypes with length <3 observed, ensure they are are 0-padded (e.g., \"008\")")
    #Add the groups
    d=addPhewasGroups(d)
    
    #Call phenotype plot as normal.
    if(add.phewas.descriptions) {
      d=addPhewasDescription(d,for.plots=T)
      phenotypeManhattan(d, annotate.phenotype.description=T,...)
    } else {
      phenotypeManhattan(d,...)
    }
  }