# Run from the Shiny app directory: Rscript check_scenarios.R
source(file.path("R", "model_bundle_helpers.R"))
source(file.path("R", "scenario_helpers.R"))
source_model_functions(app_dir = getwd())

# Existing bundles store baseline columns as matrices. A baseline compared with
# itself must produce zero changes without dropping geometry or reordering tracts.
for (path in list_available_bundles("data")) {
  bundle <- readRDS(path)
  result <- bundle$tract_data |>
    left_join(extract_counterfactual_table(bundle, bundle$baseline), by = "GEOID") |>
    add_comparison_columns()

  change_cols <- names(result)[grepl("^(delta|pct)_", names(result))]
  stopifnot(
    length(change_cols) == 22L,
    identical(result$GEOID, bundle$tract_data$GEOID),
    identical(st_geometry(result), st_geometry(bundle$tract_data)),
    all(vapply(result[change_cols] |> st_drop_geometry(), function(x) {
      is.numeric(x) && is.null(dim(x)) && all(is.na(x) | abs(x) < 1e-12)
    }, logical(1)))
  )
  cat("Baseline comparisons passed:", basename(path), "\n")
}

# Keep undefined percentage changes missing, and preserve the existing
# fractional-change convention (0.1 means 10%). Check vector and matrix inputs.
for (matrix_input in c(FALSE, TRUE)) {
  example <- result[1:5, ]
  baseline <- c(10, 0, NA_real_, -10, 20)
  counterfactual <- c(11, 5, 2, -8, NA_real_)
  example$bl_w <- if (matrix_input) matrix(baseline, ncol = 1) else baseline
  example$cf_w <- if (matrix_input) matrix(counterfactual, ncol = 1) else counterfactual
  example <- add_comparison_columns(example)
  stopifnot(
    isTRUE(all.equal(example$delta_w, c(1, 5, NA_real_, 2, NA_real_))),
    isTRUE(all.equal(example$pct_w, c(0.1, NA_real_, NA_real_, -0.2, NA_real_))),
    is.null(dim(example$pct_w))
  )
}

# Exercise a real solve, output joins, and summary generation with a tract shock.
bundle <- load_region_bundle("data", 7)
interventions <- build_intervention_record(
  tract_ids = bundle$model_inputs$GEOID[[1]],
  tract_set_name = "regression_check",
  target_variable = "b",
  shock_percent = 10,
  order = 1
)
scenario <- run_counterfactual(bundle, interventions)
stopifnot(
  nrow(scenario$result_sf) == bundle$model_inputs$N,
  nrow(scenario$summary_table) == 6L,
  all(is.finite(scenario$summary_table$counterfactual)),
  all(is.finite(scenario$result_sf$pct_w)),
  any(abs(scenario$result_sf$delta_w) > 1e-10)
)
cat("Scenario regression checks passed.\n")
