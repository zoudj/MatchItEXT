#' Compute the ratio of the variances of propensity score in the two groups
#'
#'This function accepts a MatchIt object (i.e., the result of matchit function)
#' , and calculates variance of propensity score in the two groups before and
#' after matching. The variance ratio is an indicator proposed by Rubin (2001)
#' to assess the similarity of distributions between groups.
#'
#' @param mi_obj A matchit object derived from MatchIt pacakge
#' @keywords variance ratio
#' @seealso
#' @return Return a vector of variances and variance ratios
#' @export
#' @examples
#' > \code{m_out <- matchit(treat ~ re74 + re75 + age + educ + hisp + black,
#' data = lalonde, method = "full")}
#' > \code{compute_var_ratio(m_out)}
#' @references Rubin, D. B. (2001). Using propensity scores to help design
#' observational studies: Application to the tobacco litigation. \emph{Health
#' Services and Outcomes Research Methodology, 2}(3/4), 169-188.
#' https://doi.org/10.1023/A:1020363010465

compute_var_ratio <- function(mi_obj){
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
#' @keywords grouping covariate
#' @return Return a list including the grouping variable and covariate vector
#' @export
#' @import dplyr
#' @examples
#' > \code{m_out <- matchit(treat ~ re74 + re75 + age + educ + hisp + black,
#' data = lalonde, method = "full")}
#' > \code{parse_formula(m_out)}
parse_formula <- function(mi_obj = NULL){
  grouping <- as.character(mi_obj$formula)[2]
  covariates_vec <- as.character(mi_obj$formula)[3] %>%
    strsplit(split = "+", fixed = TRUE) %>%
    unlist() %>%
    trimws(which = "both")
  return_list <- list(grouping = grouping, covariates_vec = covariates_vec)
  return(return_list)
}
