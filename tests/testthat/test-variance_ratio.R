lal <- MatchIt::lalonde

formu <- as.formula(treat ~ re74 + re75 + age + educ + black + hispan)
m_near <- MatchIt::matchit(formula = formu, data = lal, distance = "logit",
                           method = "nearest")

test_that("compute var ratio", {
  dat <- data.frame()
  expect_error(compute_var_ratio(dat),
               "The input data is not a matchit object.")
})

test_that("compute res var ratio", {
  expect_error(compute_res_var_ratio(lal, m_near,
    type_vec = c(1,2,3)), message(strwrap("The length of type vector is not
    equal to the length of covariate vector.")))
  expect_warning(compute_res_var_ratio(lal, m_near,
    type_vec = c(4, 1, 1, 1, 2, 2)), message(strwrap("There's something wrong
    in your sepecification of covariate types at vector element 1.")))
  expect_warning(compute_res_var_ratio(lal, m_near,
    type_vec = c(0, 1, 5, 1, 2, 2)), message(strwrap("There's something wrong
    in your sepecification of covariate types at vector element 3.")))
})

