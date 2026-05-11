library(tidyverse)
library(sf)

source_model_functions <- function() {
  source(file.path("R", "model_helper_functions.R"))
  source(file.path("R", "model_solver_functions.R"))
}

load_region_bundle <- function(data_dir, region_id) {
  readRDS(file.path(data_dir, paste0("region_", region_id, "_bundle.rds")))
}

list_available_bundles <- function(data_dir) {
  list.files(
    path = data_dir,
    pattern = "^region_[0-9]+_bundle\\.rds$",
    full.names = TRUE
  ) |>
    sort()
}
