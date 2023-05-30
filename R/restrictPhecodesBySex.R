restrictPhecodesBySex <- function(phenotypes,id.sex, sex.restriction=PheWAS::sex_restriction) {
  data=merge(phenotypes,id.sex,by=1,all.x=T)
  #Get the column of the sex
  g=dim(data)[2]
  #Get the restrictions found in the phenotypes data frame
  current_sex_restriction=sex.restriction[sex.restriction$phecode %in% names(phenotypes)[-1],]
  #Get male and female-only phenotypes
  male_only=current_sex_restriction[current_sex_restriction$male_only,"phecode"]
  female_only=current_sex_restriction[current_sex_restriction$female_only,"phecode"]
  #Set row column matches to NA where inds of a gender meet restricted phenotypes
  data[!is.na(data[,g])&data[,g]!="F",female_only]=NA
  data[!is.na(data[,g])&data[,g]!="M",male_only]=NA

  #Return everything, sans sex
  data[,-g]
}
