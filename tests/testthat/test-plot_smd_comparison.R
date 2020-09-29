test_that("plot smd", {
  dat <- data.frame()
  expect_error(plot_smd(dat),
    "This is not a data set derieved from MatchItEXT::compute_smd().")
})

test_that("plot ps qq",{
  dat <- list()
  expect_error(plot_ps_qq(dat), "The input is not a matchit object!")
})
