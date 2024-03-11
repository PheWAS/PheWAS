#' Add phecode descriptions and group information to existing data.
#'
#' \code{addPhecodeInfo} adds the phecode description and group information,
#' found in the data frame \code{pheinfo}. Requires the full three digit
#' phecodes in character format. This function replaces addPhewasDescription
#'  and addPhewasGroups.
#'
#'This function provides a quick way to add phecode descriptions and group
#'names.
#' It will select the first column matching any of "phecode", "phewas" or
#' "pheno" to merge on. If none is identified, it reports a warning and attempts
#'  the merge with the first column.
#' Phecodes must be character vectors, otherwise they may have inaccurate
#' mappings. If PheWAS codes are factors it will convert them to characters,
#' ≥≤give a warning, and attempt to map them.
#' If a character (or factor) variable or vector is passed in, it will return
#' code descriptions for those codes.
#'
#' @param data Data frame containing the a column with phecodes. Can
#' additionally be a character vector of phecodes. See details for requirements.
#' @param descriptions Add the phecode descriptions? Default is TRUE (yes).
#' @param groups Add the phecode group names? Default is TRUE (yes).
#' @param groupnums dd the phecode group numbers? Default is FALSE (no). Used in
#'  plotting
#' @param groupcolors Add the phecode group colors? Default is FALSE (no). Used
#' in plotting
#'
#' @return Data frame with added columns:
#' \item{description}{The description of the phecode.}
#' \item{group}{The name of the phecode group.}
#' \item{groupnum}{The assigned number of the phecode group.}
#' \item{color}{The default plotting color of the phecode group.}
#' @export
#' @import PheWASmaps
#' @importFrom utils methods
#' @importFrom methods is
#'
#' @examples 

#' addPhecodeInfo(example_phewas)
addPhecodeInfo <- function(data, descriptions=T, groups=T, groupnums=F, groupcolors=F) {
  #Convert a vector of phecodes to a data frame
 # if(class(data)[1] %in% c("character", "factor")) {data=data.frame(phenotype=data,stringsAsFactors=F)}
  #only accept data frames.
  if(class(data)[1] != 'data.frame'){stop('Data is not a Data Frame')}
  names=names(data)
 # print(names)
 # print('hi')
  

  #Find the likely phecode column
  first_match=grep("phenotype",names,ignore.case=T)[1]
  #print('hi')
  if(is.na(first_match)) {
    stop("Name matching 'pheno' not found.")
   # name=names[1]
  } else {
    name=names[first_match]
  }
#print(name)
  #Check to make sure the selected column is the correct class
  if (!is(data[[name]], 'character')) {
    
      stop("Non-character phenotypes passed in, so an accurate phecode mapping is not possible.")
    
  }

  data = inner_join(data, PheWASmaps::pheinfo, by = setNames("phecode", name))

  if(!descriptions) data = data %>% select(-description)
  if(!groups) data = data %>% select(-group)
  if(!groupnums) data = data %>% select(-groupnum)
  if(!groupcolors) data = data %>% select(-color)

  data
}
