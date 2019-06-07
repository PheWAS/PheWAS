addPhecodeInfo <- function(data, descriptions=T, groups=T, groupnums=F, groupcolors=F) {
  #Convert a vector of phecodes to a data frame
  if(class(data)[1] %in% c("character", "factor")) {data=data.frame(phenotype=data,stringsAsFactors=F)}
  names=names(data)
  #Find the likely phecode column
  first_match=grep("pheno|phewas|phecode",names,ignore.case=T)[1]
  if(is.na(first_match)) {
    warning("Name matching 'pheno', 'phecode', or 'phewas' not found, using the first column")
    name=names[1]
  } else {
    name=names[first_match]
  }
  
  #Check to make sure the selected column is the correct class
  if (class(data[[name]]) != "character") {
    if (class(data[[name]]) == "factor") {
      warning("Factor phenotype input mapped to characters")
      data[,name]=as.character(data[,name])
    } else {
      stop("Non-character or non-factor phenotypes passed in, so an accurate phecode mapping is not possible.")
    }
  }

  data = inner_join(data, pheinfo, by = setNames("phecode", name))
  
  if(!descriptions) data = data %>% select(-description)
  if(!groups) data = data %>% select(-group)
  if(!groupnums) data = data %>% select(-groupnum)
  if(!groupcolors) data = data %>% select(-color)
  
  data
}
