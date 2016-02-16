
dist_pearson <- function(x, y) {
  suppressWarnings({
    1 - cor(x, y, method = 'pearson')
  })
}

dist_euclid <- function(x, y) {
  sqrt(sum((x - y) ^ 2))
}

#' Make all possible templates for short time series gene expression data.
#' @param timepoints An integer between 3 and 9.
#' @param magnitude An integer between 1 and 4.
#' @return A matrix with one row per time point and one column per template.
#' @export
make_templates <- function(timepoints = 3, magnitude = 1) {
  # seq(-1, 1, by = 0.5)
  # [1] -1.0 -0.5  0.0  0.5  1.0
  if (!is.wholenumber(timepoints) || timepoints > 9 || timepoints < 3) {
    stop("timepoints must be a whole number between 3 and 9")
  }
  if (!is.wholenumber(magnitude) || magnitude > 4 || magnitude < 1) {
    stop("magnitude must be a whole number between 1 and 4")
  }
  magnitudes <- seq(-magnitude, magnitude)
  all_templates <- gtools::permutations(
    n = length(magnitudes),
    r = timepoints - 1,
    v = magnitudes,
    repeats.allowed = TRUE
  )
  all_templates <- cbind(rep(0, nrow(all_templates)), all_templates)
  all_templates <- apply(all_templates, 1, cumsum)
  colnames(all_templates) <- 1:ncol(all_templates)
  all_templates
}

#' Make all possible templates for short time series gene expression data.
#' @param all_templates A matrix, as returned by \code[custard]{make_templates}
#' @param templates An integer.
#' @return A matrix with one row per time point and one column per template.
#' @export
distinct_templates <- function(all_templates, templates = 10) {

  all_templates[, round(seq(1, ncol(all_templates), length.out = templates))]

  # # Include the first profile.
  # selected <- c("1")
  # selected_templates <- list(all_templates[, selected])
  #
  # # Remove it from the pool of available profiles.
  # not_selected <- !colnames(all_templates) %in% selected
  # all_templates <- all_templates[, not_selected]
  #
  # # Randomize the order of the available profiles.
  # all_templates <- all_templates[,sample(1:ncol(all_templates))]
  #
  # # Add n-1 more profiles to our selection.
  # for (i in 2:templates) {
  #   distances <- apply(all_templates, 2, function(y) {
  #     min(sapply(selected_templates, function(p) {
  #       # dist_pearson(y, p)
  #       dist_euclid(y, p)
  #     }))
  #   })
  #   next_profile <- names(distances)[which.max(distances)]
  #
  #   selected <- c(selected, next_profile)
  #   selected_templates[[i]] <- all_templates[, next_profile]
  #
  #   not_selected <- !colnames(all_templates) %in% selected
  #   all_templates <- all_templates[, not_selected]
  # }
  # do.call(cbind, selected_templates)
}
