#' Compute standardized mean differences before and after matching
#'
#' This function accepts a MatchIt object (i.e., the result of matchit function)
#' , and calculates standardized mean differences before and after matching.
#' Note exact matching and subclassification are not applicable to this function
#' . In addition, SMD can be calculated on the basis of the standard deviation
#' of original treatment group, which is the formula used in matchit function,
#' or on the basis of the simple pooled standard deviation of original treatment
#' and control group. The default is sd = "pooled", but it can be switched to
#' "treatment".
#'
#' @param mi_obj A matchit object derived from MatchIt pacakge
#' @param sd The standard deviation used as the denominator in the formula
#' @keywords SMD
#' @aliases SMD
#' @return Return a data frame containing SMD and other information
#' @export
#' @examples
#' calculate SMD based on the simple pooled standard deviation
#' > \code{smd_data <- compute_smd(mi_obj, sd = "pooled")}
#' calculate SMD based on standard deviation of original treatment group
#' > \code{smd_data <- compute_smd(mi_obj, sd = "treatment")}
#' take lalonde data as an example
#' run \code{\link[MatchIt]{matchit}} to obtain the matching result (i.e., a
#' matchit object)
#' > \code{m_out <- matchit(treat ~ re74 + re75 + age + educ + hisp + black,
#' data = lalonde, method = "optimal")}
#' use matching result and \code{\link{compute_smd}} to obtain a SMD data
#' frame
#' > \code{opt_smd <- compute_smd(m_out, sd = "treatment")}
#'
#' @references Austin, P. C. (2011). An Introduction to Propensity Score Methods
#'   for Reducing the Effects of Confounding in Observational Studies.
#'   \emph{Multivariate Behavioral Research, 46}(3), 399-424.
#'   https://doi.org/10.1080/00273171.2011.568786
#' @references Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
#'   Nonparametric Preprocessing for Parametric Causal Inference. \emph{Journal
#'   of Statistical Software, 42}(8). https://doi.org/10.18637/jss.v042.i08
#'
compute_smd<- function(mi_obj, sd = "pooled"){
  # check matchit object
  # check matching method, result of exact matching or subclassification is not
  # applicable to this function
  if (is(mi_obj, "matchit") == FALSE) {
    warning("The input data is not a matchit object.")
  } else if (as.character(mi_obj$call)[4] == "exact" |
             as.character(mi_obj$call)[4] == "subclass") {
      warning("The matching method is either exact matching or subclassification
              , SMD table cannot be generated.")
  }
    else {
      summary_data <- summary(mi_obj)
      covariate_data <- as.data.frame(mi_obj$X)
      smd_data <- data.frame(
        "mean_tr_bf" = summary_data$sum.all$`Means Treated`,
        "mean_ctl_bf" = summary_data$sum.all$`Means Control`,
        "mean_tr_af" = summary_data$sum.matched$`Means Treated`,
        "mean_ctl_af" = summary_data$sum.matched$`Means Control`)
      smd_data$rnames <- row.names(summary_data$sum.all)
      row.names(smd_data) <- row.names(summary_data$sum.all)
      smd_data$var_type[smd_data$rnames == "distance"] <- "continuous"
      for (i in 2: nrow(smd_data)){
        if (dim(table(covariate_data[row.names(smd_data)[i]])) == 2){
          smd_data$var_type[i] <- "binary"
        } else {
          smd_data$var_type[i] <- "continuous"
        }
      }
      smd_data$var_tr_bf <- NA
      smd_data$var_ctl_bf <- NA
      tr_idx <- as.numeric(as.character(mi_obj$treat))
      smd_data$var_tr_bf[1] <- var(mi_obj$distance[tr_idx == 1])
      smd_data$var_ctl_bf[1] <- var(mi_obj$distance[tr_idx == 0])
      for (i in 2:nrow(smd_data)) {
        if (smd_data$var_type[i] == "continuous") {
          smd_data$var_tr_bf[i] <-
            var(covariate_data[row.names(smd_data)[i]][tr_idx == 1, ])
          smd_data$var_ctl_bf[i] <-
            var(covariate_data[row.names(smd_data)[i]][tr_idx == 0, ])
        } else if (smd_data$var_type[i] == "binary") {
          smd_data$var_tr_bf[i] <-
            smd_data$mean_tr_bf[i] * (1 - smd_data$mean_tr_bf[i])
          smd_data$var_ctl_bf[i] <-
            smd_data$mean_ctl_bf[i] * (1 - smd_data$mean_ctl_bf[i])
        }
      }
      smd_data$sd_bf_pooled <-
        sqrt((smd_data$var_tr_bf + smd_data$var_ctl_bf)/2)
      smd_data$sd_bf_tr <- sqrt(smd_data$var_tr_bf)
      if (sd == "pooled") {
        smd_data$smd_before <- (smd_data$mean_tr_bf - smd_data$mean_ctl_bf) /
          smd_data$sd_bf_pooled
        smd_data$smd_after <-  (smd_data$mean_tr_af - smd_data$mean_ctl_af) /
          smd_data$sd_bf_pooled
      } else if (sd == "treatment") {
        smd_data$smd_before <- (smd_data$mean_tr_bf - smd_data$mean_ctl_bf) /
          smd_data$sd_bf_tr
        smd_data$smd_after <-  (smd_data$mean_tr_af - smd_data$mean_ctl_af) /
          smd_data$sd_bf_tr
      }
    }
  smd_data$rnames <- NULL
  return(smd_data)
}




