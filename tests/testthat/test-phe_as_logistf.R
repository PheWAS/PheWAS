test_that("Function works", {
  expect_equal(lapply(full_list, FUN = PheWAS:::phe_as_logistf, addidive.genotypes = TRUE, conf.level = NA, my.data = test_data_phenotype_2, min.records = 20, return.moodels = F, factor.contrasts = PheWAS:::contr.phewas, strata = NA), phe_as_log_test_1)
})
