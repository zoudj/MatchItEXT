#' Draw a dot-line plot that compares the SMD change before and after matching
#'
#' This function accepts the data set from compute_smd() and draws a dot-line
#' plot to compare the SMD change before and after matching. This function
#' differs from plot(summary(data, standardize = T) in MatchIt package in the
#' following aspects: (1) it uses the original SMD values, rather than the
#' absolute SMD values; (2) overall SMD line is colored in blue; (3) increased
#' SMD line is colored in brick red; (4) decreased SMD line is colored in gray.
#' This function depends on ggplot2. If users are not satisfied with the plot,
#' it can be revised with ggplot2. Relevant data and codes are stored in the
#' returned list.
#'
#' @param smd_data a data frame derived from \code{\link{compute_smd}}
#' @keywords SMD plot
#' @return Return a list of relevant data, code, and plot
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' > \code{plot_smd_data <- plot_smd(smd_data)}
plot_smd <- function(smd_data = NULL) {
  if (class(smd_data)[2] != "smd.data") {
    warning("This is not a data set derieved from MatchItEXT::compute_smd(). ")
  } else if (class(smd_data)[2] == "smd.data"){
      class(smd_data) <- class(smd_data)[1]
      df <- smd_data %>%
        dplyr::select(smd_before, smd_after) %>%
        tibble::rownames_to_column("var") %>%
        dplyr::mutate(trend =
                        ifelse(abs(smd_after) < abs(smd_before), "decrease",
                               ifelse(abs(smd_after) > abs(smd_before),
                                      "increase", "no change")))
      df_increase <- df %>%
        dplyr::filter(trend == "increase")
      df <- df %>%
        tidyr::pivot_longer(!c(var, trend), names_to = "stage",
                            values_to = "smd")

      df$trend[df$var == "distance"] <- "overall"
      df$stage <- factor(df$stage, levels = c("smd_before", "smd_after"))
      df$trend <- factor(df$trend, levels =
                           c("overall", "decrease", "increase", "no change"))
      colors <- c("royalblue", "gray70", "firebrick", "springgreen")

      p_smd<- ggplot(data = df, aes(x = stage, y = smd, group = var)) +
        labs(title = "Comparison of SMD Before and After Matching",
             x = NULL, y = "Standardized Mean Difference" ) +
        geom_hline(yintercept = c(-0.2, -0.1, 0, 0.1, 0.2), color = "gray90") +
        scale_y_continuous(breaks = round(seq(-0.3, 1, 0.1), digits = 1)) +
        scale_x_discrete(breaks=c("smd_before","smd_after"),
                         labels=c("All Data", "Matched Data")) +
        geom_line(aes(colour = trend)) +
        geom_point(aes(colour = trend)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_colour_manual(values=colors)
      p_smd_code <- quote(ggplot(data = df,
              aes(x = stage, y = smd, group = var)) +
              labs(title = "Comparison of SMD Before and After Matching",
                            x = NULL, y = "Standardized Mean Difference" ) +
              geom_hline(yintercept = c(-0.2, -0.1, 0, 0.1, 0.2),
                         color = "gray90") +
              scale_y_continuous(breaks = round(seq(-0.3, 1, 0.1),
                                                digits = 1)) +
              scale_x_discrete(breaks=c("smd_before","smd_after"),
                                        labels=c("All Data", "Matched Data")) +
              geom_line(aes(colour = trend)) +
              geom_point(aes(colour = trend)) +
              theme_bw() +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              scale_colour_manual(values=colors))
      results<- list(data = df, smd_increase = df_increase,
                     plot_code = p_smd_code, colors = colors, plot = p_smd)
      return(results)
  }

}


