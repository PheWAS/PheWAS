#data <- sample_data
 #phenotype_data <- createPhenotypes(data$id.vocab.code.count, id.sex = data$id.sex)
#final_data <- dplyr::inner_join(dplyr::inner_join(data$id.sex, data$genotypes), 
#phenotype_data)
#test_phewas <- phewas_ext(names(phenotype_data)[-1], 
#genotypes = c('rsEXAMPLE'), covariates = 'sex', data = final_data)

#test_that('Function Works',{
#expect_equal(addPhecodeInfo(test_phewas), donut)})  
#tests to do 
#Weird inputs
#not data frame
#test_phewas2 <- as.vector(test_phewas)
library(PheWASmaps)
test_that('Non Data Frame Produces Error', {
expect_error(addPhecodeInfo(test_phewas_no_dataframe), 'Data is not a Data Frame')})
#test_phewas3 <- as.data.table(test_phewas)
#colnames(test_phewas3)[colnames(test_phewas3) == 'phenotype'] <- 'donut'
test_that('Non phenotype column Produces Error', {
  expect_error(addPhecodeInfo(as.data.frame(test_phewas_no_pheno)), "Name matching 'pheno' not found.")
})
#Phenotype not character
test_that('Non Character phenotype column Produces Error', {
  expect_error(addPhecodeInfo(test_phewas_no_chara), 'Non-character phenotypes passed in, so an accurate phecode mapping is not possible.')
})
#Class checks
#output checks
#description
test_that('-Description works', {
  expect_equal(addPhecodeInfo(test_phewas, descriptions = F), testPhecodeInfoDescriptions)})
#group 
test_that('-Group works', {
  expect_equal(addPhecodeInfo(test_phewas, groups = F), testPhecodeInfoGroup)})

test_that('-Groupnum works', {
  expect_equal(addPhecodeInfo(test_phewas, groupnums = F), testPhecodeInfoGroupnum)})

test_that('-Groupcolor works', {
  expect_equal(addPhecodeInfo(test_phewas, groupcolors  = F), testPhecodeInfoGroupcolor)})