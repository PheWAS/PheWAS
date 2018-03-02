restrictPhecodesByGender <- function(phenotypes,id.gender) {
  data=merge(phenotypes,id.gender,by=1,all.x=T)
  #Get the column of the gender
  g=dim(data)[2]
  #Get the restrictions found in the phenotypes data frame
  current_gender_restriction=gender_restriction[gender_restriction$phecode %in% names(phenotypes)[-1],]
  #Get male and female-only phenotypes
  male_only=current_gender_restriction[current_gender_restriction$male_only,"phecode"]
  female_only=current_gender_restriction[current_gender_restriction$female_only,"phecode"]
  #Set row column matches to NA where inds of a gender meet restricted phenotypes
  data[!is.na(data[,g])&data[,g]!="F",female_only]=NA
  data[!is.na(data[,g])&data[,g]!="M",male_only]=NA

  #Return everything, sans gender
  data[,-g]
}
