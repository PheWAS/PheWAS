# test_that("Working Test Phewas", {
#   expect_equal(phenotypeManhattan(addPhecodeInfo(test_phewas, groupnums = T, groupcolors = T)), test_phenotype_manhattan_1)
#   
# })
# test_that("Missing P", {
#   expect_error(phenotypeManhattan(addPhecodeInfo(test_phewas_no_p, groupnums = T, groupcolors = T)), "Data input must contain columns phenotype and p.")
# })