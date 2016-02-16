#' Test enrichment of profiles in short time course gene expression data.
#' @param mat A matrix with gene expression profiles (columns).
#' @param templates An integer, or else a matrix of templates.
#' @param magnitude At each time point, a gene may change up to this many units.
#'   Only used if \code{templates} is an integer.
#' @param alpha Significance level (uncorrected for multiple testing).
#'   Default is 0.05.
#' @param permutations Compute this many permutations of the time points.
#' @param mc.cores Use this many cores to compute the permutations.
#' @return A list with results.
#' @export
custard <- function(
  mat,
  templates = 50,
  magnitude = 1,
  alpha = 0.05,
  permutations = 1000,
  mc.cores = 4
) {

  retval <- list()

  timepoints <- nrow(mat)

  if (timepoints > 9) {
    stop("CUSTARD is not designed to handle more than 9 time points")
  } else if (timepoints < 3) {
    stop("CUSTARD is not designed to handle fewer than 3 time points")
  }

  if (length(templates) == 1 && is.wholenumber(templates)) {
    all_templates <- make_templates(timepoints, magnitude)
    selected_templates <- distinct_templates(all_templates, templates)
    # Scale the templates to match the range of 95% of the data.
    mat_range <- abs(
      Reduce("-", quantile(mat[nrow(mat),], probs = c(0.025, 0.975)))
    )
    selected_templates <- selected_templates / mat_range
  } else if (is.matrix(templates) || is.data.frame(templates)) {
    selected_templates <- templates
  } else {
    stop("'templates' must be an integer or a matrix of templates (columns)")
  }

  template_sds <- apply(selected_templates, 2, sd)
  selected_templates <- selected_templates[, template_sds > 0]

  # Copy time point labels from the gene expression profiles.
  rownames(selected_templates) <- rownames(mat)

  retval[["templates"]] <- as.data.frame(selected_templates)

  # Compute 1 - Pearson between each gene and each template.
  cors <- 1 - cor(mat, selected_templates)

  # Assign each gene to one or more templates.
  matches <- apply(cors, 1, function(gene) {
    x <- gene == min(gene)
    x / sum(x)
  })
  observed <- rowSums(matches)
  retval[["matches"]] <- matches

  # Shuffle time points and check how many genes are assigned to each template.
  perms <- gtools::permutations(timepoints, timepoints)
  if (nrow(perms) < 1e3) {
    nulls <- apply(perms, 1, function(idx) {
      cors <- 1 - cor(mat[idx, ], selected_templates)
      matches <- apply(cors, 1, function(gene) {
        x <- gene == min(gene)
        x / sum(x)
      })
      rowSums(matches)
    })
  } else {
    perms_list <- split(perms, 1:nrow(perms))[1:1e3]
    nulls <- parallel::mclapply(perms_list, function(idx) {
      cors <- 1 - cor(mat[idx, ], selected_templates)
      matches <- apply(cors, 1, function(gene) {
        x <- gene == min(gene)
        x / sum(x)
      })
      rowSums(matches)
    }, mc.cores = mc.cores)
    nulls <- do.call(cbind, nulls)
  }
  retval[["nulls"]] <- nulls

  # Expected number of genes assigned to each template.
  expected <- rowSums(nulls) / ncol(nulls)

  # Binomial P values.
  binom_pvals <- pbinom(
    q = round(observed),
    size = ncol(mat),
    prob = expected / ncol(mat),
    lower.tail = FALSE
  )

  # Exact permutation P values.
  perm_pvals <- sapply(1:nrow(nulls), function(i) {
    (sum(nulls[i,] >= observed[i]) + 1) / (ncol(nulls) + 1)
  })

  # Bonferroni alpha level corrected for multiple testing.
  bonf_alpha <- 1 - (1 - alpha) ^ (1 / length(expected))

  # Binomial quantiles to help visualize which results are significant.
  binom_qvals <- data.frame(
    "Expected" = seq(0, max(expected), length.out = 100)
  )
  binom_qvals[["BinomQ"]] <- qbinom(
    p = 1 - bonf_alpha,
    size = ncol(mat),
    prob = binom_qvals[["Expected"]] / ncol(mat)
  )
  retval[["binom_qvals"]] <- binom_qvals

  # Binomial quantiles that correspond to each expected value.
  # binom_qvals <- qbinom(
  #   p = 1 - bonf_alpha,
  #   size = ncol(mat),
  #   prob = expected / ncol(mat)
  # )

  # Dataframe of results describing each template.
  results <- data.frame(
    "Expected" = expected,
    "Observed" = observed,
    "BinomP" = binom_pvals,
    # "BinomSignif" = observed > binom_qvals,
    "BinomSignif" = binom_pvals < bonf_alpha,
    "PermP" = perm_pvals,
    "PermSignif" = perm_pvals < bonf_alpha
  )
  retval[["results"]] <- results[order(results$BinomP), ]

  return(retval)
}
