#' A function to spread a long table of phecodes to a wide table
#' 
#' The function takes an input of 3 columns - person_ID, phecode, and count.
#' It returns a wide table with phecodes as TRUE/FALSE/NA. It can
#'  optionally use the PheWAS exclusion criteria.  The function is meant to be 
#'  used as a helper function for createPhenotype, but it can be used on its own
#'  as well
#'
#' @param phemapped a data table of 3 columns: person_ID, phecode, count
#' @param min.code.count The minimum code count to be considered a case. NA
#' results in a continuous output.
#' @param add.phecode.exclusions Apply PheWAS exclusions to phecodes.
#' @param id.sex If supplied, restrict the phecodes by sex. This should be a
#' data frame with the first column being the id and the second the sex, "M" or
#' "F", of the individual. Individuals with any other specification will have
#' all sex specific phenotypes set to NA.
#' @param full.population.ids List of IDs in the "complete" population. This
#' allows for individuals with no observed codes to have appropriate "control"
#' status, eg 0s or FALSE in every field
#' @param aggregate.fun Aggregate function for duplicated phenotypes
#' (phecodes, etc) in an individual. The default supports a naive
#' "distinct date" approach. Use \code{sum} to support count data
#' @param vocabulary.map Map between supplied vocabularies and phecodes. Allows
#' for custom phecode maps. By default uses
#' \code{\link[PheWASmaps:phecode_map]{PheWASmaps::phecode_map}},
#' which supports ICD9CM (v1.2) and ICD10CM (beta-2018).
#'  The package also includes the ICD10 beta map
#'  (\code{\link[PheWASmaps:phecode_map_icd10]{PheWASmaps::phecode_map_icd10}}), which
#'  can be used in this parameter.
#' @param rollup.map Map between phecodes and all codes that they expand to, eg
#'  parent codes. By default uses the PheWASmaps::phecode_rollup_map.
#' @param exclusion.map Map between phecodes and their exclusions. By default
#'  uses the PheWASmaps::phecode_exclude.
#' @param gender.exclusion Map determining sex-based exclusions
#' @param id.name the name of the phecodes
#'
#' @return A data frame. The first column contains the supplied id for each
#'  individual (preserving the name of the original column). The following
#'   columns are all present phewas codes. They contain T/F/NA for
#'   case/control/exclude or continuous/NA if min.code.count was NA.
#'   
#' @export
#'
#' @examples
#' \donttest{
#' pheSpread(PheWAS:::test_phemapped, full.population.ids=
#' unique(sample_data$id.vocab.code.index[[1]]), id.name = 
#' names(sample_data$id.vocab.code.count)[1], aggregate.fun = sum, 
#' id.sex = sample_data$id.sex)}
pheSpread <- function(phemapped, min.code.count=2,add.phecode.exclusions = T, 
                      id.sex,
                      full.population.ids=unique(id.vocab.code.index[[1]]),
                      aggregate.fun=default_code_agg,
                      vocabulary.map=PheWASmaps::phecode_map,
                      rollup.map=PheWASmaps::phecode_rollup_map,
                      exclusion.map=PheWASmaps::phecode_exclude,
                      gender.exclusion = PheWASmaps::gender_restriction,
                      id.name){
  message("Aggregating codes...")
  phecode=ungroup(summarize(group_by(phemapped,id,code),count=aggregate.fun(index)))
  phecode=phecode[phecode$count>0,]
  #Check exclusions, and add them to the list
  if(add.phecode.exclusions) {
    message("Mapping exclusions...")
    exclusions = inner_join(phecode %>% rename(exclusion_criteria=code), exclusion.map, by = "exclusion_criteria")
    exclusions = exclusions %>%  transmute(id, code, count=-1) %>% distinct()
    phecode=rbind(phecode,exclusions)
  }
  #If there is request for a min code count, adjust counts to -1 if needed
  if(!is.na(min.code.count)&(max(!is.na(phecode$count)&phecode$count<min.code.count))) {
    phecode[!is.na(phecode$count)&phecode$count<min.code.count,]$count=-1
  }
  if(!is.na(min.code.count)|add.phecode.exclusions) {
    message("Coalescing exclusions and min.code.count as applicable...")
    phecode=ungroup(summarize(group_by(phecode,id,code),count=max(count)))
  }
  message("Reshaping data...")
  phens=spread(phecode,code,count,fill=0)
  #Set exclusions to NA, preserving IDs just in case one is -1
  tmp_id=phens[,1]
  phens[phens==-1]=NA
  phens[,1]=tmp_id
  #Add in inds present in input or the full population list, but without mapped phecodes
  missing_ids=setdiff(full.population.ids,phens[["id"]])
  if(length(missing_ids)>0) {
    empty_record=phens[1,-1]
    empty_record[]=0
    phens=rbind(phens,data.frame(id=missing_ids,empty_record,check.names=F))
  }
  #Change to logical if there is a min code count
  if(!is.na(min.code.count)) {phens[,-1]=phens[,-1]>0}
  #If there are sex restrictions, set them to NA
  if(!missing(id.sex)) {
    phens=restrictPhecodesBySex(phens,id.sex, gender.exclusion)
  }
  #Limit to full population ids
  phens = filter(phens, id %in% full.population.ids)
  #Rename the ID column to the input ID column name
  if(!missing(id.name)){
  names(phens)[1]=id.name
  }
  phens
}