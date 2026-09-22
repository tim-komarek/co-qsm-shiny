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
Land supply in those bundles should be built from `K_available_tracts`, not raw `ALAND`.
If `K_available_tracts` is exactly zero, the bundle builder keeps the displayed value at zero and applies a tiny positive floor only inside the solver inputs.

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

Run `Rscript check_scenarios.R` from the app directory to check cached-bundle
comparisons and a counterfactual solve before deploying.

- `app.R`: Shiny app entry point
- `app.css`: responsive layout for desktop embeds and phones
- `app.js`: map resizing and collapsible phone controls
- `R/`: model and app helper code
- `data/`: precomputed region bundles used by the app

## Website embedding and Connect Cloud

The desktop app fills the height and width of its window or iframe. Model
controls scroll independently, and the Run scenario button stays visible.
The tract table has its own tab; scenario results have Map, Summary, and
Interventions views. Tables scroll within their panels instead of widening
the page. Below 768 pixels, controls collapse and the page scrolls vertically.
There is no REDI website header in the app.

Deploy this folder with `app.R` as the primary file. Both `app.css` and `app.js`
are included in `manifest.json` and must be published along with the R files
and regional bundles. No additional R packages are required for the layout.

Use the public app URL in the website iframe:
`https://019e1840-3311-c8d9-c976-d1ba24f5de17.share.connect.posit.cloud/`.
Set the iframe to `width: 100%; height: 100%; border: 0; display: block` and give
its containing element an explicit height. The website's current 700-pixel
minimum height can still exceed a short browser window; adjust that rule in
WordPress if the whole embed needs to fit below the site's header and title.
The Shiny app cannot change the parent website's header, footer, or iframe size.

## Notes

- This is a standalone Shiny-app repository, not the full `CO_QSM` project.
- The repo includes the minimal model code needed to run the app locally.
