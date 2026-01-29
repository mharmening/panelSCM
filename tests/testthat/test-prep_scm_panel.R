test_that("prep_scm_panel returns correct structure", {
  # Skip on CRAN to avoid long running tests
  skip_on_cran()

  # Create simple test data
  df <- data.frame(
    unit_id = rep(1:10, each = 10),
    unit_name = rep(paste0("Unit_", 1:10), each = 10),
    time = rep(1:10, 10),
    outcome = rnorm(100),
    pred1 = rnorm(100),
    treatment = 0
  )
  df$treatment[df$unit_id <= 3 & df$time >= 7] <- 1

  result <- prep_scm_panel(
    df = df,
    outcome_var = "outcome",
    predictor_vars = "pred1",
    id_var = "unit_id",
    name_var = "unit_name",
    time_var = "time",
    treatment_indicator = "treatment",
    donor_range = 3,
    parallel = FALSE
  )

  expect_type(result, "list")
  expect_true("optimal_n" %in% names(result))
  expect_true("models" %in% names(result))
})

