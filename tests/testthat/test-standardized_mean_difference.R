lal <- MatchIt::lalonde

formu <- as.formula(treat ~ re74 + re75 + age + educ + black + hispan)
m_exact <- MatchIt::matchit(formula = formu, data = lal, method = "exact")
m_sub <- MatchIt::matchit(formula = formu, data = lal, distance = "logit",
                          method = "subclass", subclass = 6)
m_near <- MatchIt::matchit(formula = formu, data = lal, distance = "logit",
                           method = "nearest")

test_that("compute smd", {
  dat <- data.frame()
  expect_error(compute_smd(dat), "The input data is not a matchit object.")
  expect_error(compute_smd(m_exact), message(strwrap("The matching method is
    either exact matching or subclassification, SMD table cannot be generated.",
                                             )))
  expect_error(compute_smd(m_sub), message(strwrap("The matching method is
    either exact matching or subclassification, SMD table cannot be generated.",
                                           )))
  expect_equal(class(compute_smd(m_near)), c("data.frame", "smd.data"))
  expect_error(compute_smd(m_near, sd= 'something'),
    "The argument of sd can only be either 'pooled' or 'treatment'.")
})

test_that("compute sub smd", {
  dat <- data.frame()
  expect_error(compute_sub_smd(dat), "The input data is not a matchit object.")
  expect_error(compute_sub_smd(m_exact), message(strwrap("The matching method
    is exact matching, compute_sub_smd() and compute_smd() are inapplicable to
    it.")))
  expect_error(compute_sub_smd(m_near), message(strwrap("The matching method is
    not subclassification, please try compute_smd() instead.")))
  expect_error(compute_sub_smd(m_sub, sd= 'something'),
    "The argument of sd can only be either 'pooled' or 'treatment'.")
})
