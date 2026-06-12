# Core counterfactual mediation methods for the proteome-COPD manuscript.
#
# Component-protein models are evaluated for proteins in the fixed signature
# after the corresponding signature-level mediation result satisfies the
# robust-positive rule.

mc_require_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_packages) > 0L) {
    stop(
      sprintf("Install required R package(s): %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

mc_check_columns <- function(data, columns, context = "input data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("%s is missing column(s): %s", context, paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

z_standardize <- function(x, center = NULL, scale = NULL) {
  x <- as.numeric(x)
  if (is.null(center)) {
    center <- mean(x, na.rm = TRUE)
  }
  if (is.null(scale)) {
    scale <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) {
    stop("Cannot z-standardize a variable with missing or non-positive scale.", call. = FALSE)
  }
  out <- (x - center) / scale
  attr(out, "center") <- center
  attr(out, "scale") <- scale
  out
}

build_mediation_data <- function(data,
                                 time_col,
                                 event_col,
                                 exposure_col,
                                 mediator,
                                 covariates = character()) {
  mediator_values <- if (is.character(mediator) && length(mediator) == 1L && mediator %in% names(data)) {
    data[[mediator]]
  } else {
    mediator
  }
  required_columns <- unique(c(time_col, event_col, exposure_col, covariates))
  mc_check_columns(data, required_columns, "mediation data")
  if (length(mediator_values) != nrow(data)) {
    stop("mediator must be a column name or a vector with one value per row.", call. = FALSE)
  }

  model_df <- data.frame(
    time = as.numeric(data[[time_col]]),
    event = as.numeric(data[[event_col]]),
    exposure = as.numeric(data[[exposure_col]]),
    mediator = as.numeric(mediator_values),
    data[, covariates, drop = FALSE],
    check.names = FALSE
  )
  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  model_df <- model_df[is.finite(model_df$time) & model_df$time > 0, , drop = FALSE]
  list(data = model_df, covariates = covariates)
}

fit_cmaverse_cox_mediation <- function(model_data, covariates = character()) {
  mc_require_packages("CMAverse")
  mc_check_columns(model_data, c("time", "event", "exposure", "mediator", covariates), "CMAverse model data")
  suppressWarnings(CMAverse::cmest(
    data = model_data,
    model = "rb",
    outcome = "time",
    event = "event",
    exposure = "exposure",
    mediator = "mediator",
    basec = covariates,
    EMint = FALSE,
    yreg = "coxph",
    mreg = list("linear"),
    astar = 0,
    a = 1,
    mval = list(0),
    estimation = "paramfunc",
    inference = "delta",
    yrare = TRUE
  ))
}

tidy_cmaverse_effects <- function(fit_obj) {
  effect_map <- c(
    Rte = "te",
    Rpnde = "nde",
    Rtnde = "tnde",
    Rpnie = "pnie",
    Rtnie = "tnie",
    Rcde = "cde",
    pm = "pm"
  )
  available <- intersect(names(effect_map), names(fit_obj$effect.pe))
  out <- data.frame(
    effect = unname(effect_map[available]),
    cmaverse_effect = available,
    estimate = as.numeric(fit_obj$effect.pe[available]),
    std_error = as.numeric(fit_obj$effect.se[available]),
    lower = as.numeric(fit_obj$effect.ci.low[available]),
    upper = as.numeric(fit_obj$effect.ci.high[available]),
    p_value = as.numeric(fit_obj$effect.pval[available]),
    stringsAsFactors = FALSE
  )
  ratio_effect <- out$effect != "pm"
  out$hr <- NA_real_
  out$log_estimate <- NA_real_
  out$log_lower <- NA_real_
  out$log_upper <- NA_real_
  out$hr[ratio_effect] <- out$estimate[ratio_effect]
  out$log_estimate[ratio_effect] <- log(out$estimate[ratio_effect])
  out$log_lower[ratio_effect] <- log(out$lower[ratio_effect])
  out$log_upper[ratio_effect] <- log(out$upper[ratio_effect])
  out
}

empty_mediation_effect_columns <- function() {
  effects <- c("te", "nde", "tnde", "pnie", "tnie", "cde", "pm")
  fields <- c("estimate", "std_error", "lower", "upper", "p_value", "hr", "log_estimate", "log_lower", "log_upper")
  values <- rep(NA_real_, length(effects) * length(fields))
  names(values) <- as.vector(outer(effects, fields, paste, sep = "__"))
  as.data.frame(as.list(values), check.names = FALSE)
}

effects_to_wide <- function(effect_df) {
  if (nrow(effect_df) == 0L) {
    return(empty_mediation_effect_columns())
  }
  fields <- c("estimate", "std_error", "lower", "upper", "p_value", "hr", "log_estimate", "log_lower", "log_upper")
  out <- list()
  for (i in seq_len(nrow(effect_df))) {
    effect_name <- effect_df$effect[[i]]
    for (field in fields) {
      out[[paste(effect_name, field, sep = "__")]] <- effect_df[[field]][[i]]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

run_signature_mediation <- function(data,
                                    time_col,
                                    event_col,
                                    exposure_col,
                                    signature_score_col,
                                    covariates = character()) {
  run_one_mediator(
    data = data,
    time_col = time_col,
    event_col = event_col,
    exposure_col = exposure_col,
    mediator = signature_score_col,
    mediator_name = signature_score_col,
    mediator_type = "signature_score",
    covariates = covariates
  )
}

run_one_mediator <- function(data,
                             time_col,
                             event_col,
                             exposure_col,
                             mediator,
                             mediator_name = "mediator",
                             mediator_type = "signature_or_component",
                             covariates = character()) {
  base_row <- data.frame(
    exposure = exposure_col,
    mediator = mediator_name,
    mediator_type = mediator_type,
    n_total = NA_integer_,
    n_events = NA_integer_,
    exposure_mean = NA_real_,
    exposure_sd = NA_real_,
    mediator_mean = NA_real_,
    mediator_sd = NA_real_,
    n_covariates = length(covariates),
    mediation_package = "CMAverse",
    mediation_model = "rb_paramfunc_delta_coxph",
    interaction_used = FALSE,
    model_status = "failed",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  built <- tryCatch(
    build_mediation_data(
      data = data,
      time_col = time_col,
      event_col = event_col,
      exposure_col = exposure_col,
      mediator = mediator,
      covariates = covariates
    ),
    error = function(e) e
  )
  if (inherits(built, "error")) {
    base_row$error_message <- conditionMessage(built)
    return(cbind(base_row, empty_mediation_effect_columns()))
  }

  model_df <- built$data
  base_row$n_total <- nrow(model_df)
  base_row$n_events <- sum(model_df$event == 1, na.rm = TRUE)
  base_row$exposure_mean <- mean(model_df$exposure, na.rm = TRUE)
  base_row$exposure_sd <- stats::sd(model_df$exposure, na.rm = TRUE)
  base_row$mediator_mean <- mean(model_df$mediator, na.rm = TRUE)
  base_row$mediator_sd <- stats::sd(model_df$mediator, na.rm = TRUE)
  if (nrow(model_df) == 0L || base_row$n_events == 0L || length(unique(model_df$event)) < 2L) {
    base_row$error_message <- "No evaluable rows or no outcome events after complete-case filtering."
    return(cbind(base_row, empty_mediation_effect_columns()))
  }

  fit <- tryCatch(
    fit_cmaverse_cox_mediation(model_data = model_df, covariates = built$covariates),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    base_row$error_message <- conditionMessage(fit)
    return(cbind(base_row, empty_mediation_effect_columns()))
  }

  effects <- tidy_cmaverse_effects(fit)
  base_row$model_status <- "ok"
  base_row$error_message <- NA_character_
  cbind(base_row, effects_to_wide(effects), stringsAsFactors = FALSE)
}

add_bh_fdr <- function(results,
                       p_col = "tnie__p_value",
                       q_col = "q_tnie",
                       group_cols = NULL) {
  out <- as.data.frame(results, stringsAsFactors = FALSE)
  if (!p_col %in% names(out)) {
    out[[q_col]] <- NA_real_
    return(out)
  }
  out[[q_col]] <- NA_real_
  if (is.null(group_cols) || length(group_cols) == 0L) {
    ok <- !is.na(out[[p_col]])
    out[[q_col]][ok] <- stats::p.adjust(out[[p_col]][ok], method = "BH")
    return(out)
  }
  mc_check_columns(out, group_cols, "FDR grouping data")
  groups <- split(seq_len(nrow(out)), interaction(out[, group_cols, drop = FALSE], drop = TRUE))
  for (idx in groups) {
    ok <- !is.na(out[[p_col]][idx])
    if (any(ok)) {
      out[[q_col]][idx[ok]] <- stats::p.adjust(out[[p_col]][idx][ok], method = "BH")
    }
  }
  out
}

add_robust_positive_flags <- function(results,
                                      require_q = FALSE,
                                      q_col = "q_tnie",
                                      q_threshold = 0.05) {
  out <- as.data.frame(results, stringsAsFactors = FALSE)
  for (column in c("pm__estimate", "pm__lower", "tnie__lower")) {
    if (!column %in% names(out)) {
      out[[column]] <- NA_real_
    }
  }
  status_ok <- if ("model_status" %in% names(out)) out$model_status == "ok" else rep(TRUE, nrow(out))
  out$pm_positive <- !is.na(out$pm__estimate) & out$pm__estimate > 0
  out$pm_ci_positive <- !is.na(out$pm__lower) & out$pm__lower > 0
  out$tnie_ci_positive <- !is.na(out$tnie__lower) & out$tnie__lower > 1
  out$robust_positive_mediation <- status_ok & out$pm_positive & out$pm_ci_positive & out$tnie_ci_positive
  if (isTRUE(require_q)) {
    if (!q_col %in% names(out)) {
      out[[q_col]] <- NA_real_
    }
    out$retain_for_interpretation <- out$robust_positive_mediation & !is.na(out[[q_col]]) & out[[q_col]] < q_threshold
  } else {
    out$retain_for_interpretation <- out$robust_positive_mediation
  }
  out
}

run_component_mediation_models <- function(data,
                                           time_col,
                                           event_col,
                                           exposure_col,
                                           mediator_cols,
                                           covariates = character(),
                                           signature_level_positive = TRUE,
                                           standardize_mediators = TRUE,
                                           adjust_fdr = TRUE,
                                           require_q = TRUE) {
  if (!isTRUE(signature_level_positive)) {
    stop(
      "Run component-protein mediation after confirming that the corresponding signature-level mediation result satisfies the robust-positive rule.",
      call. = FALSE
    )
  }
  mc_check_columns(data, mediator_cols, "component mediation data")
  results <- do.call(rbind, lapply(mediator_cols, function(mediator_col) {
    mediator_values <- if (isTRUE(standardize_mediators)) {
      z_standardize(data[[mediator_col]])
    } else {
      as.numeric(data[[mediator_col]])
    }
    run_one_mediator(
      data = data,
      time_col = time_col,
      event_col = event_col,
      exposure_col = exposure_col,
      mediator = mediator_values,
      mediator_name = mediator_col,
      mediator_type = if (isTRUE(standardize_mediators)) "z_standardized_component_protein" else "component_protein",
      covariates = covariates
    )
  }))
  if (isTRUE(adjust_fdr)) {
    results <- add_bh_fdr(results)
  }
  add_robust_positive_flags(results, require_q = require_q)
}
