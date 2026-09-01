library(tidyverse)
library(sf)

source_model_functions <- function(app_dir = getwd()) {
  source(file.path(app_dir, "R", "model_helper_functions.R"))
  source(file.path(app_dir, "R", "model_solver_functions.R"))
}

bundle_file_path <- function(data_dir, region_id) {
  region_id <- as.integer(region_id)
  file.path(data_dir, paste0("region_", region_id, "_bundle.rds"))
}

load_region_bundle <- function(data_dir, region_id) {
  readRDS(bundle_file_path(data_dir = data_dir, region_id = region_id))
}

list_available_bundles <- function(data_dir) {
  list.files(
    path = data_dir,
    pattern = "^region_[0-9]+_bundle\\.rds$",
    full.names = TRUE
  ) |>
    sort()
}

read_tract_data <- function(project_dir) {
  data_path <- file.path(project_dir, "1_CO_data", "df_CO_QSM.RData")

  if (!file.exists(data_path)) {
    stop("Could not find `1_CO_data/df_CO_QSM.RData`.")
  }

  tract_df <- readRDS(data_path)

  required_cols <- c(
    "GEOID",
    "CountyName",
    "region",
    "region_name",
    "res",
    "emp",
    "ALAND",
    "K_available_tracts",
    "Q_index",
    "geometry"
  )

  missing_cols <- setdiff(required_cols, names(tract_df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The tract file is missing: ",
        paste(missing_cols, collapse = ", "),
        "."
      )
    )
  }

  sf::st_as_sf(tract_df)
}

build_region_bundle <- function(project_dir, region_id, app_dir = file.path(project_dir, "Shiny")) {
  region_id <- as.integer(region_id)

  if (is.na(region_id)) {
    stop("`region_id` must be coercible to an integer.")
  }

  source_model_functions(app_dir = app_dir)

  tract_sf <- read_tract_data(project_dir)
  region_sf <- tract_sf[tract_sf$region == region_id, ]

  if (nrow(region_sf) == 0) {
    stop(paste0("Could not find region ", region_id, " in the tract file."))
  }

  osrm_path <- file.path(project_dir, "1_CO_data", "orsm_tracts_by_region.rds")

  if (!file.exists(osrm_path)) {
    stop("Could not find `1_CO_data/orsm_tracts_by_region.rds`.")
  }

  osrm_tracts <- readRDS(osrm_path)
  region_key <- paste0("region_", region_id)

  if (!region_key %in% names(osrm_tracts)) {
    stop(paste0("Could not find travel times for ", region_key, "."))
  }

  time_df <- osrm_tracts[[region_key]]$time_df
  row_order <- match(time_df$location_id, region_sf$GEOID)

  if (any(is.na(row_order))) {
    stop("Travel-time rows could not be matched to tract GEOIDs.")
  }

  region_sf <- region_sf[row_order, ]
  time_cols <- paste0("time_", region_sf$GEOID)

  if (!all(time_cols %in% names(time_df))) {
    stop("Travel-time columns do not line up with tract GEOIDs.")
  }

  t_ij <- as.matrix(time_df[, time_cols])
  N <- nrow(region_sf)
  L_j <- region_sf$emp
  L_i <- region_sf$res
  L_i <- L_i * sum(L_j) / sum(L_i)
  K_raw <- region_sf$K_available_tracts
  # The updated developable-land measure can be exactly zero in a few tracts.
  # Keep those zeros in the displayed data, but use a tiny positive floor in
  # the solver inputs so the fixed-point iteration remains numerically defined.
  K <- pmax(K_raw, 1)
  Q <- region_sf$Q_index

  inversion <- inversionModel(
    N = N,
    L_i = L_i,
    L_j = L_j,
    Q = Q,
    K = K,
    t_ij = t_ij
  )

  baseline <- solveModel(
    N = N,
    L_i = L_i,
    L_j = L_j,
    K = K,
    t_ij = t_ij,
    a = inversion$a,
    b = inversion$b,
    varphi = inversion$varphi,
    w_eq = inversion$w,
    u_eq = inversion$u,
    Q_eq = inversion$Q_norm,
    ttheta_eq = inversion$ttheta
  )

  tract_data <- region_sf |>
    mutate(
      res_obs = res,
      emp_obs = emp,
      Q_obs = Q_index,
      K_obs = K_available_tracts,
      a_bl = inversion$a,
      b_bl = inversion$b,
      varphi_bl = inversion$varphi,
      inv_A = inversion$A,
      inv_B = inversion$B,
      inv_w = inversion$w,
      inv_u = inversion$u,
      inv_Q_norm = inversion$Q_norm,
      inv_ttheta = inversion$ttheta,
      bl_w = baseline$w,
      bl_W_i = baseline$W_i,
      bl_B = baseline$B,
      bl_A = baseline$A,
      bl_Q = baseline$Q,
      bl_L_i = baseline$L_i,
      bl_L_j = baseline$L_j,
      bl_ybar = baseline$ybar,
      bl_lambda_i = baseline$lambda_i,
      bl_ttheta = baseline$ttheta,
      bl_u = baseline$u
    )

  list(
    metadata = tibble(
      region_id = region_id,
      region_name = region_sf$region_name[[1]],
      n_tracts = N,
      baseline_U = inversion$U[[1]],
      baseline_solve_U = baseline$U[[1]]
    ),
    tract_data = tract_data,
    model_inputs = list(
      N = N,
      GEOID = region_sf$GEOID,
      L_i = L_i,
      L_j = L_j,
      Q = Q,
      K = K,
      t_ij = t_ij
    ),
    inversion = list(
      a = inversion$a,
      b = inversion$b,
      varphi = inversion$varphi,
      w = inversion$w,
      u = inversion$u,
      Q_norm = inversion$Q_norm,
      ttheta = inversion$ttheta,
      A = inversion$A,
      B = inversion$B,
      U = inversion$U
    ),
    baseline = list(
      w = baseline$w,
      W_i = baseline$W_i,
      B = baseline$B,
      A = baseline$A,
      Q = baseline$Q,
      L_i = baseline$L_i,
      L_j = baseline$L_j,
      ybar = baseline$ybar,
      lambda_i = baseline$lambda_i,
      ttheta = baseline$ttheta,
      u = baseline$u,
      U = baseline$U
    )
  )
}

write_region_bundle <- function(bundle, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_file <- bundle_file_path(
    data_dir = output_dir,
    region_id = bundle$metadata$region_id[[1]]
  )

  saveRDS(bundle, out_file)
  out_file
}
