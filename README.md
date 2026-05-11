# Colorado Quantitative Spatial Model Shiny App

This repository contains a standalone Shiny app for exploring tract-level counterfactuals from the Colorado quantitative spatial model.

## What the app does

- loads precomputed baseline model bundles for 8 Colorado regions
- displays baseline tract maps
- lets users select arbitrary tract sets
- applies stacked tract-level percent shocks
- solves a new equilibrium and maps counterfactual changes

The app currently supports tract-level percent shocks to:

- `a`: productivity
- `b`: amenities
- `varphi`: development density
- `K`: land supply
- `Q`: floorspace prices
- `w`: wages
- `L_i`: residents
- `L_j`: workers

## Included data

This repo includes precomputed regional `.rds` bundles in `data/`, so the app should run immediately after clone.

## Required R packages

The app expects these packages to be installed:

- `shiny`
- `leaflet`
- `tidyverse`
- `sf`

## Run locally

Open the repository in R or RStudio and run:

```r
shiny::runApp()
```

## Repository structure

- `app.R`: Shiny app entry point
- `R/`: model and app helper code
- `data/`: precomputed region bundles used by the app

## Notes

- This is a standalone Shiny-app repository, not the full `CO_QSM` project.
- The repo includes the minimal model code needed to run the app locally.
