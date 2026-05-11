library(tidyverse)
library(sf)

policy_variable_choices <- function() {
  c(
    "Productivity (a)" = "a",
    "Amenity (b)" = "b",
    "Development Density (varphi)" = "varphi",
    "Land Supply (K)" = "K"
  )
}

baseline_variable_choices <- function() {
  c(
    "Residents (observed)" = "res_obs",
    "Employment (observed)" = "emp_obs",
    "Floorspace Price (observed)" = "Q_obs",
    "Land Supply (observed)" = "K_obs",
    "Recovered Productivity (a)" = "a_bl",
    "Recovered Amenity (b)" = "b_bl",
    "Recovered Development Density (varphi)" = "varphi_bl",
    "Baseline Wage" = "bl_w",
    "Baseline Floorspace Price" = "bl_Q",
    "Baseline Residents" = "bl_L_i",
    "Baseline Employment" = "bl_L_j",
    "Baseline Welfare" = "bl_u"
  )
}

result_variable_choices <- function() {
  c(
    "Wages" = "w",
    "Floorspace Prices" = "Q",
    "Residents" = "L_i",
    "Employment" = "L_j",
    "Welfare" = "u",
    "Productivity (effective A)" = "A",
    "Amenities (effective B)" = "B",
    "Commercial Floorspace Share" = "ttheta"
  )
}

display_mode_choices <- function() {
  c(
    "Counterfactual Level" = "cf",
    "Absolute Change" = "delta",
    "Percent Change" = "pct"
  )
}

empty_interventions <- function() {
  tibble(
    order = integer(),
    tract_set_name = character(),
    tract_ids = list(),
    target_variable = character(),
    shock_percent = numeric(),
    tract_count = integer()
  )
}

build_intervention_record <- function(tract_ids, tract_set_name, target_variable, shock_percent, order) {
  tibble(
    order = order,
    tract_set_name = tract_set_name,
    tract_ids = list(sort(unique(tract_ids))),
    target_variable = target_variable,
    shock_percent = shock_percent,
    tract_count = length(unique(tract_ids))
  )
}

apply_interventions <- function(bundle, interventions) {
  scenario_inputs <- list(
    a = bundle$inversion$a,
    b = bundle$inversion$b,
    varphi = bundle$inversion$varphi,
    K = bundle$model_inputs$K
  )

  if (nrow(interventions) == 0) {
    return(scenario_inputs)
  }

  tract_lookup <- tibble(
    GEOID = bundle$model_inputs$GEOID,
    row_id = seq_along(bundle$model_inputs$GEOID)
  )

  interventions |>
    arrange(order) |>
    pwalk(function(order, tract_set_name, tract_ids, target_variable, shock_percent, tract_count) {
      row_ids <- tract_lookup |>
        filter(GEOID %in% tract_ids) |>
        pull(row_id)

      multiplier <- 1 + (shock_percent / 100)
      scenario_inputs[[target_variable]][row_ids] <<- scenario_inputs[[target_variable]][row_ids] * multiplier
    })

  scenario_inputs
}

extract_counterfactual_table <- function(bundle, solve_output) {
  tibble(
    GEOID = bundle$model_inputs$GEOID,
    cf_w = as.numeric(solve_output$w),
    cf_W_i = as.numeric(solve_output$W_i),
    cf_B = as.numeric(solve_output$B),
    cf_A = as.numeric(solve_output$A),
    cf_Q = as.numeric(solve_output$Q),
    cf_L_i = as.numeric(solve_output$L_i),
    cf_L_j = as.numeric(solve_output$L_j),
    cf_ybar = as.numeric(solve_output$ybar),
    cf_lambda_i = as.numeric(solve_output$lambda_i),
    cf_ttheta = as.numeric(solve_output$ttheta),
    cf_u = as.numeric(solve_output$u)
  )
}

add_comparison_columns <- function(result_sf) {
  comparison_pairs <- tribble(
    ~baseline_col, ~counterfactual_col, ~base_name,
    "bl_w", "cf_w", "w",
    "bl_W_i", "cf_W_i", "W_i",
    "bl_B", "cf_B", "B",
    "bl_A", "cf_A", "A",
    "bl_Q", "cf_Q", "Q",
    "bl_L_i", "cf_L_i", "L_i",
    "bl_L_j", "cf_L_j", "L_j",
    "bl_ybar", "cf_ybar", "ybar",
    "bl_lambda_i", "cf_lambda_i", "lambda_i",
    "bl_ttheta", "cf_ttheta", "ttheta",
    "bl_u", "cf_u", "u"
  )

  output <- result_sf

  for (i in seq_len(nrow(comparison_pairs))) {
    baseline_col <- comparison_pairs$baseline_col[[i]]
    counterfactual_col <- comparison_pairs$counterfactual_col[[i]]
    base_name <- comparison_pairs$base_name[[i]]

    output[[paste0("delta_", base_name)]] <- output[[counterfactual_col]] - output[[baseline_col]]
    output[[paste0("pct_", base_name)]] <- dplyr::if_else(
      abs(output[[baseline_col]]) > 0,
      (output[[counterfactual_col]] - output[[baseline_col]]) / output[[baseline_col]],
      NA_real_
    )
  }

  output
}

build_summary_table <- function(bundle, result_sf, solve_output) {
  tibble(
    metric = c(
      "Aggregate welfare (U)",
      "Mean tract welfare (u)",
      "Total residents (L_i)",
      "Total employment (L_j)",
      "Mean wage (w)",
      "Mean floorspace price (Q)"
    ),
    baseline = c(
      bundle$baseline$U,
      mean(result_sf$bl_u, na.rm = TRUE),
      sum(result_sf$bl_L_i, na.rm = TRUE),
      sum(result_sf$bl_L_j, na.rm = TRUE),
      mean(result_sf$bl_w, na.rm = TRUE),
      mean(result_sf$bl_Q, na.rm = TRUE)
    ),
    counterfactual = c(
      as.numeric(solve_output$U),
      mean(result_sf$cf_u, na.rm = TRUE),
      sum(result_sf$cf_L_i, na.rm = TRUE),
      sum(result_sf$cf_L_j, na.rm = TRUE),
      mean(result_sf$cf_w, na.rm = TRUE),
      mean(result_sf$cf_Q, na.rm = TRUE)
    )
  ) |>
    mutate(
      absolute_change = counterfactual - baseline,
      percent_change = if_else(abs(baseline) > 0, absolute_change / baseline, NA_real_)
    )
}

run_counterfactual <- function(bundle, interventions) {
  scenario_inputs <- apply_interventions(bundle, interventions)

  solve_output <- solveModel(
    N = bundle$model_inputs$N,
    L_i = bundle$model_inputs$L_i,
    L_j = bundle$model_inputs$L_j,
    K = scenario_inputs$K,
    t_ij = bundle$model_inputs$t_ij,
    a = scenario_inputs$a,
    b = scenario_inputs$b,
    varphi = scenario_inputs$varphi,
    w_eq = bundle$inversion$w,
    u_eq = bundle$inversion$u,
    Q_eq = bundle$inversion$Q_norm,
    ttheta_eq = bundle$inversion$ttheta
  )

  result_sf <- bundle$tract_data |>
    left_join(extract_counterfactual_table(bundle, solve_output), by = "GEOID") |>
    add_comparison_columns()

  list(
    result_sf = result_sf,
    solve_output = solve_output,
    summary_table = build_summary_table(bundle, result_sf, solve_output)
  )
}

resolve_result_column <- function(base_name, display_mode) {
  if (display_mode == "cf") {
    return(paste0("cf_", base_name))
  }

  if (display_mode == "delta") {
    return(paste0("delta_", base_name))
  }

  paste0("pct_", base_name)
}
