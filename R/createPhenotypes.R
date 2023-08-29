#' Creates a phenotype table from id, ICD9CM, ICD10CM (or phecode, etc), data.
#'
#' This function takes a data frame with four columns: id, vocabulary_id, code,
#'  and index. It returns a wide table with phecodes as TRUE/FALSE/NA. It can
#'  optionally use the PheWAS exclusion criteria.
#'
#' By default, this function returns a wide format data frame with boolean
#' phenotypes suitable for PheWAS analysis. Specifying a \code{min.code.count=NA}
#' will permit continuous code count phenotypes.
#' The default exclusions can be skipped with \code{add.exclusions=F}. In
#'  conjuntion with \code{translate=F} (and optionally adjusting
#'  \code{min.code.count} and \code{aggregate.fun}), one can use this function
#'  as a simple reshaping wrapper.
#' @param id.vocab.code.index Data frame with four columns of information: id,
#' vocabulary_id, code, and index. The id and index columns can have other
#' names, but id must be consistent among input files. The vocabulary_id and
#' code must match up with the vocabulary.map file. The default supports the
#' vocabularies "ICD9CM" and "ICD10CM". Code contains the raw code value. Note
#' the beta ICD10 map is provided in
#' \code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}},
#'  which can be provided as a parameter to
#' \code{vocabulary.map}.
#' @param min.code.count The minimum code count to be considered a case. NA
#' results in a continuous output.
#' @param add.phecode.exclusions Apply PheWAS exclusions to phecodes.
#' @param translate Should the input be translated to phecodes? Defaults to
#' TRUE. Generally recommended, though can be skipped if phecodes are provided.
#' @param id.sex If supplied, restrict the phecodes by sex. This should be a
#' data frame with the first column being the id and the second the sex, "M" or
#' "F", of the individual. Individuals with any other specification will have
#' all sex specific phenotypes set to NA.
#' @param full.population.ids List of IDs in the "complete" population. This
#' allows for individuals with no observed codes to have appropriate "control"
#' status, eg 0s or FALSE in every field.
#' @param aggregate.fun Aggregate function for duplicated phenotypes
#' (phecodes, etc) in an individual. The default supports a naive
#' "distinct date" approach. Use \code{sum} to support count data.
#' @param vocabulary.map Map between supplied vocabularies and phecodes. Allows
#' for custom phecode maps. By default uses
#' \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}},
#' which supports ICD9CM (v1.2) and ICD10CM (beta-2018).
#'  The package also includes the ICD10 beta map
#'  (\code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}}), which
#'  can be used in this parameter.
#' @param rollup.map Map between phecodes and all codes that they expand to, eg
#'  parent codes. By default uses the PheWAS::phecode_rollup_map.
#' @param exclusion.map Map between phecodes and their exclusions. By default
#'  uses the PheWAS::phecode_exclude.
#'
#' @return A data frame. The first column contains the supplied id for each
#'  individual (preserving the name of the original column). The following
#'   columns are all present phewas codes. They contain T/F/NA for
#'   case/control/exclude or continuous/NA if min.code.count was NA.
#'
#' @export
#'
#' @examples id_icd9_count=data.frame(id=c(1,2,2,2,3),vocabulary_id="icd9",
#' code=c("714","250.11","714.1","714","250.11"),
#' count=c(1,5,1,1,0))
#' createPhenotypes(id_icd9_count)
#' \donttest{
#' #Complex example
#' ex=generateExample(n=500,hit="335")
#' #Extract the two relevant parts from the returned list
#' id.vocab.code.count=ex$id.vocab.code.count
#' id.sex=ex$id.sex
#' #Create the phecode table- translates the codes, adds
#' #exclusions, and reshapes to a wide format.
#' #Sum up the counts in the data where applicable.
#' phenotypes=createPhenotypes(id.vocab.code.count,
#'                             aggregate.fun=sum, id.sex=id.sex)
#' 
#' }
createPhenotypes <-
  function(id.vocab.code.index, min.code.count=2, add.phecode.exclusions=T, translate=T, id.sex,
           full.population.ids=unique(id.vocab.code.index[[1]]),
           aggregate.fun=PheWAS:::default_code_agg,
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           exclusion.map=PheWAS::phecode_exclude)
  {
    id.name=names(id.vocab.code.index)[1]

    #Warn if id.sex information is not provided.
    if(missing(id.sex)) { warning("It is recommended to provide id.sex information to help address spurious sex-specific associations.") }

    if(!translate) {
      #Warn about exclusions if input is not translated and not phecodes. Same with id.sex
      if(add.phecode.exclusions & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")}
      if(!missing(id.sex) & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")}
      phemapped=tbl_df(data.frame(id=id.vocab.code.index[[1]],code=id.vocab.code.index[[3]],index=id.vocab.code.index[[4]],stringsAsFactors = F))
    } else {
      #check to make sure numeric codes were not passed in
      if(!class(id.vocab.code.index[[3]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id.vocab.code.index)=c("id","vocabulary_id","code","index")
      message("Mapping codes to phecodes...")
      phemapped=mapCodesToPhecodes(id.vocab.code.index, vocabulary.map=vocabulary.map, rollup.map=rollup.map) %>% transmute(id, code=phecode, index)
    }

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
      phens=restrictPhecodesBySex(phens,id.sex)
    }

    #Limit to full population ids
    phens = filter(phens, id %in% full.population.ids)

    #Rename the ID column to the input ID column name
    names(phens)[1]=id.name

    phens
  }
