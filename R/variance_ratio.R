#' Compute the ratio of the variances of propensity score in the two groups
#'
#'This function accepts a MatchIt object (i.e., the result of matchit function)
#' , and calculates variance of propensity score in the two groups before and
#' after matching. The variance ratio is an indicator proposed by Rubin (2001)
#' to assess the similarity of distributions between groups.
#'
#' @param mi_obj A matchit object derived from MatchIt pacakge
#' @keywords variance ratio
#' @seealso compute_res_var_ratio()
#' @return Return a vector of variances and variance ratios
#' @export
#' @examples
#' m_out <- MatchIt::matchit(treat ~ re74 + re75 + age + educ + hispan +
#' black, data = MatchIt::lalonde, method = "nearest")
#' compute_var_ratio(m_out)
#' @references Rubin, D. B. (2001). Using propensity scores to help design
#' observational studies: Application to the tobacco litigation. \emph{Health
#' Services and Outcomes Research Methodology, 2}(3/4), 169-188.
#' https://doi.org/10.1023/A:1020363010465

compute_var_ratio <- function(mi_obj = NULL){
  if (is(mi_obj, "matchit") == FALSE) {
    stop("The input data is not a matchit object.")
  }
  var_tr_bf <- var(mi_obj$distance[mi_obj$treat == 1], na.rm = T)
  var_ctl_bf <- var(mi_obj$distance[mi_obj$treat == 0], na.rm = T)
  ratio_bf <- var_tr_bf / var_ctl_bf
  data_af <- data.frame(distance = mi_obj$distance, wt = mi_obj$weights, tr =
                          mi_obj$treat)
  data_af_tr <- data_af[which(data_af$wt != 0 & data_af$tr == 1),]
  data_af_ctl <- data_af[which(data_af$wt != 0 & data_af$tr == 0),]
  var_tr_af <- Hmisc::wtd.var(data_af_tr$distance, weights = data_af_tr$wt,
                              method = 'ML')
  var_ctl_af <- Hmisc::wtd.var(data_af_ctl$distance, weights = data_af_ctl$wt,
                               method = 'ML')
  ratio_af <- var_tr_af / var_ctl_af
  results <- c(var_tr_bf = var_tr_bf, var_ctl_bf = var_ctl_bf, var_tr_af =
                 var_tr_af, var_ctl_af = var_ctl_af, ratio_bf = ratio_bf,
                 ratio_af = ratio_af)
  return(results)
}

#' Parse formula to obtain grouping variable and covariate vector
#'
#' This function parses the formula used in matchit() to obtain the grouping
#' variable and the covariate vector.
#'
#' @param mi_obj A matchit object derived from MatchIt pacakge
#' @keywords grouping covariate
#' @return Return a list including the grouping variable and covariate vector
#' @export
#' @import dplyr
#' @examples
#' m_out <- MatchIt::matchit(treat ~ re74 + re75 + age + educ + hispan +
#' black, data = MatchIt::lalonde, method = "nearest")
#' parse_formula(m_out)
parse_formula <- function(mi_obj = NULL){
  grouping <- as.character(mi_obj$formula)[2]
  covec <- as.character(mi_obj$formula)[3] %>%
    strsplit(split = "+", fixed = TRUE) %>%
    unlist() %>%
    trimws(which = "both")
  return_list <- list(grouping = grouping, covec = covec)
  return(return_list)
}

#' Compute ratio of variances of residuals for covariates
#'
#' This function computes the ratio of variances of residuals for covariates,
#' which was proposed by Rubin (2001). Applicable covariate types include
#' continuous, binary and ordinal. Multinomial variables are not applicable to
#' this function due to the absence of single residual. Usually a k-category
#' multinomial variable will have k-1 residuals if multinomial logistic or
#' probit regression is applied. For continuous variable, glm(family= gaussian)
#'  is used; for binary variable, glm(family= binomial(link= logit)) is used;
#' for ordinal variable, MASS::polr(method = "logistic") is used, then single
#' residual is obtained by using sure::resids().
#'
#' @keywords residual variance ratio
#' @param original_data A data frame containing original data
#' @param mi_obj A matchit object derived from MatchIt pacakge
#' @param type_vec A vector specifying covariate types, valid values: 'ordinal',
#' '3', 3; 'binary', '2', 2; 'continuous', '1', 1; 'excluded', '0', 0, NA. The
#' last one means not to compute the ratio for this covariate. The length of
#' this vector should be the same as that of covariate vector used in
#' propensity score estimation.
#' @param discard A logical value. TRUE means some observations are discarded
#' before matching (with respect to discarded argument in matchit function),
#' then the ratio before matching is based on the data after discard; FALSE
#' means no observation is discarded before matching, then the ratio before
#' matching is based on the original intact data.
#' @importFrom stats as.formula glm gaussian binomial
#' @export
#' @seealso parse_formula() compute_var_ratio()
#' @examples
#'  m_out <- MatchIt::matchit(treat ~ re74 + re75 + age + educ + hispan + black,
#'   data = MatchIt::lalonde, method = "nearest")
#' # use parse_formula() to check grouping variable and covariates
#'  parse_formula(m_out)
#'  compute_res_var_ratio(original_data = MatchIt::lalonde, mi_obj =
#'   m_out, type_vec = c(0, 1, 1, 1, 2, 2))
#' @references Rubin, D. B. (2001). Using propensity scores to help design
#' observational studies: Application to the tobacco litigation. \emph{Health
#' Services and Outcomes Research Methodology, 2}(3/4), 169-188.
#' https://doi.org/10.1023/A:1020363010465
compute_res_var_ratio<- function(original_data = NULL, mi_obj = NULL, type_vec
                                 = NULL, discard = FALSE){
  covec <- parse_formula(mi_obj)$covec
  if (length(type_vec) != length(covec)){
    stop(message(strwrap("The length of type vector is not equal to the length
      of covariate vector.")))
  } else {
    grouping_v <- as.character(mi_obj$formula)[2]

    # calculate before-matching ratio
    ratio_vec_bf <- rep(NA, length(covec))
    original_data$distance <- mi_obj$distance
    original_data$discarded <- mi_obj$discarded
    if (discard == FALSE){
      dat_bf <- original_data
    } else {
      dat_bf <- original_data[which(original_data$discarded == FALSE),]
    }
    for (k in seq_along(type_vec)){
      formu <- as.formula(paste(covec[k], "distance", sep = "~"))
      if (type_vec[k] == 0 | type_vec[k] == '0' | type_vec[k] == 'excluded' |
          is.na(type_vec[k])){
        next #excluded covariate
      } else if (type_vec[k] == 1 | type_vec[k] == '1' | type_vec[k] ==
                 'continuous'){
        # regular linear regression for continuous variable
        mod <- glm(formula = formu, family = gaussian(link = "identity"),
                   data = dat_bf)
        dat_bf$residuals <- mod$residuals
        ratio_vec_bf[k] <- with(dat_bf, var(residuals[get(grouping_v) == 1])) /
          with(dat_bf, var(residuals[get(grouping_v) == 0]))
        dat_bf$residuals <- NULL
      } else if (type_vec[k] == 2 | type_vec[k] == '2' | type_vec[k] ==
                 'binary') {
        # logistic regression for binary variable
        mod <- glm(formula = formu, family = binomial(link = "logit"), data =
                     dat_bf)
        dat_bf$residuals <- mod$residuals
        ratio_vec_bf[k] <- with(dat_bf, var(residuals[get(grouping_v) == 1])) /
          with(dat_bf, var(residuals[get(grouping_v) == 0]))
        dat_bf$residuals <- NULL
      } else if (type_vec[k] == 3 | type_vec[k] == '3' | type_vec[k] ==
                 'ordinal'){
        # ordinal logistics regression for ordinal variable
        dat_bf[,grouping_v] <- as.factor(dat_bf[,grouping_v])
        mod <- MASS::polr(formula = formu, data = dat_bf, method = "logistic")
        dat_bf$residuals <- sure::resids(mod)
        ratio_vec_bf[k] <- with(dat_bf, var(residuals[get(grouping_v) == 1])) /
          with(dat_bf, var(residuals[get(grouping_v) == 0]))
        dat_bf$residuals <- NULL
      } else {
        warning(strwrap(sprintf("There's something wrong in your sepecification
          of covariate types at vector element %d.", k)))
      }
    }

    # calculate afte-matching ratio
    ratio_vec_af <- rep(NA, length(covec))
    dat_af <- match.data(mi_obj)
    for (i in seq_along(type_vec)){
      formu <- as.formula(paste(covec[i], "distance", sep = "~"))
      if (type_vec[i] == 0 | type_vec[i] == '0' | type_vec[i] == 'excluded' |
          is.na(type_vec[i])){
        next #excluded covariate
      } else if (type_vec[i] == 1 | type_vec[i] == '1' | type_vec[i] ==
                 'continuous'){
        # regular linear regression for continuous variable
        mod <- glm(formula = formu, family = gaussian(link = "identity"),
                   data = dat_af, weights = .data$weights)
        dat_af$residuals <- mod$residuals
        ratio_vec_af[i] <- with(dat_af, var(residuals[get(grouping_v) == 1])) /
          with(dat_af, var(residuals[get(grouping_v) == 0]))
        dat_af$residuals <- NULL
      } else if (type_vec[i] == 2 | type_vec[i] == '2' | type_vec[i] ==
                 'binary'){
        # logistic regression for binary variable
        mod <- glm(formula = formu, family = binomial(link = "logit"), data =
                     dat_af, weights = .data$weights)
        dat_af$residuals <- mod$residuals
        ratio_vec_af[i] <- with(dat_af, var(residuals[get(grouping_v) == 1])) /
          with(dat_af, var(residuals[get(grouping_v) == 0]))
        dat_af$residuals <- NULL
      } else if (type_vec[i] == 3 | type_vec[i] == '3' | type_vec[i] ==
                 'ordinal'){
        # ordinal logistics regression for ordinal variable
        dat_af[,grouping_v] <- as.factor(dat_af[,grouping_v])
        mod <- MASS::polr(formula = formu, data = dat_af, weights = .data$weights,
                          method = "logistic")
        dat_af$residuals <- sure::resids(mod)
        ratio_vec_af[i] <- with(dat_af, var(residuals[get(grouping_v) == 1])) /
          with(dat_af, var(residuals[get(grouping_v) == 0]))
        dat_af$residuals <- NULL
      }
    }
  }
  return(list(before_matching = ratio_vec_bf, after_matching = ratio_vec_af))
}
