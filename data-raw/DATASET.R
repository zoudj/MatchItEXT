## code to prepare `DATASET` dataset goes here
dat <- MatchIt::lalonde

formu <- as.formula(treat ~ re74 + re75 + age + educ + black + hispan)
m_exact <- MatchIt::matchit(formula = formu, data = dat, method = "exact")
m_sub <- MatchIt::matchit(formula = formu, data = dat, distance = "logit",
                          method = "subclass", subclass = 6)
m_near <- MatchIt::matchit(formula = formu, data = dat, distance = "logit",
                           method = "nearest")

usethis::use_data(m_exact, m_sub, m_near, overwrite = TRUE)
