# Core lasso Cox protein-signature methods for the proteome-COPD manuscript.
#
# The manuscript signature model uses a lasso-penalized Cox model
# (glmnet alpha = 1), keeps age and sex unpenalized, and applies the
# lambda.min cross-validation rule.

pc_require_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_packages) > 0L) {
    stop(
      sprintf("Install required R package(s): %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pc_check_columns <- function(data, columns, context = "input data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("%s is missing column(s): %s", context, paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pc_check_complete_rows <- function(data, columns, context = "input data") {
  missing_rows <- !stats::complete.cases(data[, columns, drop = FALSE])
  if (any(missing_rows)) {
    stop(
      sprintf(
        "%s contains missing values in %d row(s). This methods code expects complete analysis inputs, matching the imputed-data setting of the manuscript analysis.",
        context,
        sum(missing_rows)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pc_numeric_matrix <- function(data, columns, context) {
  out <- as.data.frame(data[, columns, drop = FALSE])
  is_numeric <- vapply(out, is.numeric, FUN.VALUE = logical(1))
  if (!all(is_numeric)) {
    stop(
      sprintf("%s must be numeric. Non-numeric column(s): %s", context, paste(columns[!is_numeric], collapse = ", ")),
      call. = FALSE
    )
  }
  out <- as.matrix(out)
  storage.mode(out) <- "double"
  out
}

fit_lasso_cox_signature <- function(data,
                                    time_col,
                                    event_col,
                                    protein_cols,
                                    unpenalized_cols = c("age_at_recruitment", "sex"),
                                    nfolds = 10,
                                    foldid = NULL,
                                    seed = NULL,
                                    standardize = TRUE,
                                    ...) {
  pc_require_packages(c("survival", "glmnet"))
  if (length(protein_cols) == 0L) {
    stop("protein_cols must contain at least one protein column.", call. = FALSE)
  }

  required_columns <- unique(c(time_col, event_col, unpenalized_cols, protein_cols))
  pc_check_columns(data, required_columns, "signature fitting data")

  model_df <- as.data.frame(data[, required_columns, drop = FALSE])
  model_df[[time_col]] <- as.numeric(model_df[[time_col]])
  model_df[[event_col]] <- as.numeric(model_df[[event_col]])
  pc_check_complete_rows(model_df, required_columns, "signature fitting data")
  if (any(!is.finite(model_df[[time_col]]) | model_df[[time_col]] <= 0)) {
    stop("signature fitting data contains non-positive or non-finite follow-up times.", call. = FALSE)
  }
  if (any(!model_df[[event_col]] %in% c(0, 1))) {
    stop("signature fitting data event column must be coded 0/1.", call. = FALSE)
  }
  if (nrow(model_df) < 10L || length(unique(model_df[[event_col]])) < 2L) {
    stop("Signature fitting needs at least 10 complete rows and two outcome states.", call. = FALSE)
  }

  protein_matrix <- pc_numeric_matrix(model_df, protein_cols, "protein_cols")
  unpenalized_matrix <- pc_numeric_matrix(
    model_df,
    unpenalized_cols,
    "unpenalized_cols"
  )

  x <- cbind(unpenalized_matrix, protein_matrix)
  penalty_factor <- c(rep(0, ncol(unpenalized_matrix)), rep(1, ncol(protein_matrix)))
  y <- survival::Surv(model_df[[time_col]], model_df[[event_col]])

  if (is.null(foldid)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(model_df)))
  } else if (length(foldid) != nrow(model_df)) {
    stop("foldid length must match the number of complete input rows.", call. = FALSE)
  }

  cv_fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = nfolds,
    foldid = foldid,
    type.measure = "deviance",
    penalty.factor = penalty_factor,
    standardize = standardize,
    keep = FALSE,
    ...
  )

  structure(
    list(
      cv_fit = cv_fit,
      lambda_rule = "lambda.min",
      alpha = 1,
      protein_cols = protein_cols,
      unpenalized_cols = unpenalized_cols,
      model_columns = colnames(x),
      penalty_factor = penalty_factor,
      foldid = foldid,
      n_used = nrow(model_df),
      n_events = sum(model_df[[event_col]] == 1)
    ),
    class = "protein_signature_fit"
  )
}

extract_signature_coefficients <- function(signature_fit,
                                           lambda = "lambda.min",
                                           protein_cols = NULL,
                                           nonzero_only = FALSE) {
  pc_require_packages("glmnet")
  cv_fit <- if (inherits(signature_fit, "protein_signature_fit")) signature_fit$cv_fit else signature_fit
  if (!identical(lambda, "lambda.min")) {
    stop("The manuscript analysis uses lambda.min; this public helper does not expose alternate lambda rules.", call. = FALSE)
  }
  if (is.null(protein_cols) && inherits(signature_fit, "protein_signature_fit")) {
    protein_cols <- signature_fit$protein_cols
  }

  coef_matrix <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))
  out <- data.frame(
    term = rownames(coef_matrix),
    coefficient = as.numeric(coef_matrix[, 1L]),
    lambda = "lambda.min",
    stringsAsFactors = FALSE
  )
  out$term_type <- ifelse(out$term %in% protein_cols, "protein", "unpenalized_covariate")
  out$in_signature <- out$term_type == "protein" & out$coefficient != 0
  out <- out[order(out$term_type, -abs(out$coefficient), out$term), , drop = FALSE]
  rownames(out) <- NULL
  if (isTRUE(nonzero_only)) {
    out <- out[out$coefficient != 0, , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

signature_protein_weights <- function(coefficients, protein_cols = NULL) {
  if (is.null(protein_cols)) {
    protein_cols <- if ("term_type" %in% names(coefficients)) {
      coefficients$term[coefficients$term_type == "protein"]
    } else {
      coefficients$term
    }
  }

  weights <- coefficients[coefficients$term %in% protein_cols & coefficients$coefficient != 0, , drop = FALSE]
  if (nrow(weights) == 0L) {
    stop("No non-zero protein coefficients were supplied.", call. = FALSE)
  }
  weights
}

signature_raw_score <- function(data, coefficients, protein_cols = NULL) {
  weights <- signature_protein_weights(coefficients, protein_cols)
  pc_check_columns(data, weights$term, "signature scoring data")
  pc_check_complete_rows(data, weights$term, "signature scoring data")
  protein_matrix <- pc_numeric_matrix(data, weights$term, "signature scoring columns")
  list(
    score = as.numeric(protein_matrix %*% weights$coefficient),
    weights = weights
  )
}

derive_signature_scaling <- function(derivation_data,
                                     coefficients,
                                     protein_cols = NULL) {
  raw <- signature_raw_score(derivation_data, coefficients, protein_cols)
  center <- mean(raw$score)
  scale <- stats::sd(raw$score)
  if (!is.finite(scale) || scale <= 0) {
    stop("Derivation-sample signature score has missing or non-positive scale.", call. = FALSE)
  }
  list(
    center = center,
    scale = scale,
    n_derivation = length(raw$score)
  )
}

calculate_signature_score <- function(data,
                                      coefficients,
                                      protein_cols = NULL,
                                      center,
                                      scale) {
  if (missing(center) || missing(scale) || !is.finite(scale) || scale <= 0) {
    stop(
      "Provide derivation-sample center and positive scale. Estimate them once in the derivation/training sample with derive_signature_scaling(), then reuse them for validation or full-sample scoring.",
      call. = FALSE
    )
  }
  raw <- signature_raw_score(data, coefficients, protein_cols)
  out <- data.frame(
    signature_score = raw$score,
    signature_score_z = (raw$score - center) / scale
  )
  attr(out, "center") <- center
  attr(out, "scale") <- scale
  attr(out, "scoring_note") <- paste(
    "When applying a fixed signature to validation or full analysis data,",
    "use the center and scale estimated in the derivation/training sample;",
    "do not re-standardize separately within each new dataset."
  )
  out
}

pc_quote_name <- function(name) {
  sprintf("`%s`", gsub("`", "", name, fixed = TRUE))
}

fit_signature_association <- function(data,
                                      time_col,
                                      event_col,
                                      score_col = "signature_score_z",
                                      covariates = character(),
                                      robust = FALSE) {
  pc_require_packages("survival")
  required_columns <- unique(c(time_col, event_col, score_col, covariates))
  pc_check_columns(data, required_columns, "signature association data")
  model_df <- as.data.frame(data[, required_columns, drop = FALSE])
  model_df[[time_col]] <- as.numeric(model_df[[time_col]])
  model_df[[event_col]] <- as.numeric(model_df[[event_col]])
  pc_check_complete_rows(model_df, required_columns, "signature association data")
  if (any(!is.finite(model_df[[time_col]]) | model_df[[time_col]] <= 0)) {
    stop("signature association data contains non-positive or non-finite follow-up times.", call. = FALSE)
  }
  if (any(!model_df[[event_col]] %in% c(0, 1))) {
    stop("signature association data event column must be coded 0/1.", call. = FALSE)
  }

  rhs <- paste(vapply(c(score_col, covariates), pc_quote_name, character(1)), collapse = " + ")
  formula_obj <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", pc_quote_name(time_col), pc_quote_name(event_col), rhs))
  fit <- survival::coxph(formula_obj, data = model_df, robust = robust, x = TRUE, y = TRUE, model = TRUE)
  summary_fit <- summary(fit)
  coef_table <- as.data.frame(summary_fit$coefficients, stringsAsFactors = FALSE)
  coef_table$term <- rownames(coef_table)
  rownames(coef_table) <- NULL
  conf_table <- as.data.frame(summary_fit$conf.int, stringsAsFactors = FALSE)
  conf_table$term <- rownames(conf_table)
  rownames(conf_table) <- NULL
  merged <- merge(coef_table, conf_table, by = "term", suffixes = c("_coef", "_ci"))
  score_row <- merged[merged$term == score_col | merged$term == paste0("`", score_col, "`"), , drop = FALSE]
  p_col <- grep("^Pr\\(", names(score_row), value = TRUE)[1]
  hr_col <- if ("exp(coef)_ci" %in% names(score_row)) "exp(coef)_ci" else "exp(coef)"

  term_summary <- data.frame(
    term = score_col,
    hazard_ratio = if (nrow(score_row) == 1L) score_row[[hr_col]] else NA_real_,
    ci_low = if (nrow(score_row) == 1L) score_row[["lower .95"]] else NA_real_,
    ci_high = if (nrow(score_row) == 1L) score_row[["upper .95"]] else NA_real_,
    p_value = if (nrow(score_row) == 1L && length(p_col) == 1L) score_row[[p_col]] else NA_real_,
    n_used = fit$n,
    n_events = fit$nevent,
    stringsAsFactors = FALSE
  )
  list(fit = fit, term = term_summary)
}
