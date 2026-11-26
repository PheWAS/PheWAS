quick_confint_data <- function() {
  set.seed(1)
  n <- 300
  data.frame(
    id = 1:n,
    case = rbinom(n, 1, 0.4) == 1,
    level = rnorm(n),
    snp = rbinom(n, 2, 0.3)
  )
}

expect_quick_confint <- function(res, level = 0.95) {
  z <- qnorm((1 - level) / 2, lower.tail = FALSE)
  logistic <- res$type == "logistic"
  lower <- res$beta - z * res$SE
  upper <- res$beta + z * res$SE
  expect_equal(res$lower.q, ifelse(logistic, exp(lower), lower))
  expect_equal(res$upper.q, ifelse(logistic, exp(upper), upper))
}

test_that("phewas() returns quick confidence intervals", {
  d <- quick_confint_data()
  res <- suppressMessages(phewas(
    phenotypes = c("case", "level"), genotypes = "snp", data = d,
    quick.confint.level = 0.95
  ))
  expect_true(all(c("lower.q", "upper.q") %in% names(res)))
  expect_setequal(res$type, c("logistic", "linear"))
  expect_quick_confint(res)
})

test_that("phewas_ext() returns quick confidence intervals with models", {
  d <- quick_confint_data()
  res <- suppressMessages(phewas_ext(
    phenotypes = c("case", "level"), genotypes = "snp", data = d,
    quick.confint.level = 0.95, return.models = TRUE
  ))
  expect_named(res, c("results", "models"))
  expect_quick_confint(res$results)
})
