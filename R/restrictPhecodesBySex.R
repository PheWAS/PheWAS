#' Add PheWAS code descriptions to existing data.
#'
#' \code{restrictPhecodesBySex} alters a table for PheWAS with phecodes, as from
#'  \code{\link[PheWAS:createPhenotypes]{createPhenotypes}}, to exclude
#'  individuals with non-applicable sexes from certain phenotypes.
#'
#' @param phenotypes   The PheWAS table to have restrictions applied. The first
#' column should be the id.
#' @param id.sex A data frame with the first column being the id and the second
#' the gender, "M" or "F", of the individual. Individuals with any other
#'
#' specification will have all gender specific phenotypes set to NA.
#' @export
#' @return The \code{phenotypes} data frame with NA values for individuals that
#' do not match the sex for sex-specific codes.
#' @examples
#' data <- sample_data
#' restrictPhecodesBySex(data$id.vocab.code.count, data$id.sex)

restrictPhecodesBySex <- function(phenotypes,id.sex) {
  data=merge(phenotypes,id.sex,by=1,all.x=T)
  #Get the column of the sex
  g=dim(data)[2]
  #Get the restrictions found in the phenotypes data frame
  current_gender_restriction=PheWASmaps::gender_restriction[PheWASmaps::gender_restriction$phecode %in% names(phenotypes)[-1],]
  #Get male and female-only phenotypes
  male_only=current_gender_restriction[current_gender_restriction$male_only,"phecode"]
  female_only=current_gender_restriction[current_gender_restriction$female_only,"phecode"]
  #Set row column matches to NA where inds of a gender meet restricted phenotypes
  data[!is.na(data[,g])&data[,g]!="F",female_only]=NA
  data[!is.na(data[,g])&data[,g]!="M",male_only]=NA

  #Return everything, sans sex
  data[,-g]
}
