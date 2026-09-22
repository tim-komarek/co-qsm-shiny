library(shiny)
library(leaflet)
library(tidyverse)
library(sf)

app_dir <- normalizePath(getwd())
data_dir <- file.path(app_dir, "data")

source(file.path(app_dir, "R", "model_bundle_helpers.R"))
source(file.path(app_dir, "R", "scenario_helpers.R"))
source_model_functions(app_dir = app_dir)

bundle_files <- list_available_bundles(data_dir)

bundle_index <- if (length(bundle_files) > 0) {
  map_dfr(bundle_files, function(path) {
    bundle <- readRDS(path)

    bundle$metadata |>
      mutate(file_path = path)
  })
} else {
  tibble(
    region_id = integer(),
    region_name = character(),
    n_tracts = integer(),
    baseline_U = numeric(),
    baseline_solve_U = numeric(),
    file_path = character()
  )
}

baseline_choices <- baseline_variable_choices()
result_choices <- result_variable_choices()
policy_choices <- policy_variable_choices()
display_choices <- display_mode_choices()

build_palette <- function(values) {
  finite_values <- as.numeric(values)[is.finite(values)]
  if (length(finite_values) == 0) {
    finite_values <- c(0, 1)
  }

  leaflet::colorNumeric(
    palette = "viridis",
    domain = finite_values,
    na.color = "#d9d9d9"
  )
}

ui <- fluidPage(
  title = "Colorado Quantitative Spatial Model",
  tags$head(
    includeCSS(file.path(app_dir, "app.css")),
    includeScript(file.path(app_dir, "app.js"))
  ),
  div(
    class = "qsm-app",
    div(class = "qsm-title", h2("Colorado Quantitative Spatial Model")),
    if (nrow(bundle_index) == 0) {
      div(
        class = "qsm-empty",
        h4("No Shiny region bundles found."),
        p("Run `Rscript Shiny/build_region_bundles.R` from the project root, then restart the app.")
      )
    } else {
      div(
        class = "qsm-layout",
        tags$aside(
          class = "qsm-sidebar",
          `aria-label` = "Model controls",
          tags$details(
            id = "qsm-controls",
            open = "open",
            tags$summary("Model controls and tract selection"),
            div(
              class = "qsm-controls-scroll",
              selectInput(
                "region_id", "Region",
                choices = set_names(bundle_index$region_id, bundle_index$region_name),
                selected = bundle_index$region_id[[1]]
              ),
              selectInput(
                "baseline_map_var", "Baseline variable",
                choices = baseline_choices, selected = "res_obs"
              ),
              h4("Selected tracts"),
              verbatimTextOutput("selected_tract_summary"),
              selectizeInput(
                "tract_table_geoid", "Add tracts from list",
                choices = NULL, multiple = TRUE
              ),
              div(
                class = "qsm-buttons",
                actionButton("add_table_selection", "Add to selected tracts"),
                actionButton("clear_selection", "Clear selected tracts")
              ),
              h4("Intervention builder"),
              textInput("tract_set_name", "Tract set name", value = "selected_tracts"),
              selectInput("target_variable", "Policy variable", choices = policy_choices, selected = "b"),
              sliderInput("shock_percent", "Percent change", min = -95, max = 100, value = 10, step = 1, post = "%"),
              div(
                class = "qsm-buttons",
                actionButton("add_intervention", "Add intervention"),
                actionButton("remove_last_intervention", "Remove last intervention"),
                actionButton("clear_interventions", "Clear interventions")
              ),
              div(class = "qsm-status", textOutput("intervention_count"))
            )
          ),
          div(
            class = "qsm-run-area",
            actionButton("run_scenario", "Run scenario", class = "btn-primary")
          )
        ),
        tags$main(
          class = "qsm-main",
          tabsetPanel(
            id = "main_view",
            tabPanel(
              "Baseline map", value = "baseline",
              div(class = "qsm-map-panel", leafletOutput("baseline_map", height = "100%"))
            ),
            tabPanel(
              "Tract table", value = "tracts",
              div(
                class = "qsm-table-panel",
                p(class = "qsm-table-note", "First 15 tracts in the region. Select any tract using the map or the tract list in Model controls."),
                div(class = "qsm-table-scroll", tabindex = "0", `aria-label` = "Tract data", tableOutput("tract_table"))
              )
            ),
            tabPanel(
              "Scenario results", value = "results",
              tabsetPanel(
                id = "results_view", type = "pills",
                tabPanel(
                  "Map", value = "map",
                  div(
                    class = "qsm-results-panel",
                    div(
                      class = "qsm-results-controls",
                      selectInput("result_map_var", "Result variable", choices = result_choices, selected = "u"),
                      selectInput("result_display_mode", "Display mode", choices = display_choices, selected = "pct"),
                      checkboxInput("hide_selected_in_results", "Hide selected tracts", value = FALSE)
                    ),
                    div(class = "qsm-map-panel", leafletOutput("results_map", height = "100%"))
                  )
                ),
                tabPanel(
                  "Summary", value = "summary",
                  div(class = "qsm-table-panel", h4("Scenario summary"),
                      div(class = "qsm-table-scroll", tabindex = "0", `aria-label` = "Scenario summary", tableOutput("summary_table")))
                ),
                tabPanel(
                  "Interventions", value = "interventions",
                  div(class = "qsm-table-panel", h4("Intervention stack"),
                      div(class = "qsm-table-scroll", tabindex = "0", `aria-label` = "Intervention stack", tableOutput("intervention_table")))
                )
              )
            )
          )
        )
      )
    }
  )
)

server <- function(input, output, session) {
  if (nrow(bundle_index) == 0) {
    return(invisible(NULL))
  }

  current_bundle <- reactiveVal(NULL)
  bundle_cache <- reactiveValues(data = list())
  selected_ids <- reactiveVal(character())
  interventions <- reactiveVal(empty_interventions())
  scenario_output <- reactiveVal(NULL)

  observeEvent(input$region_id, {
    region_key <- as.character(input$region_id)

    if (is.null(bundle_cache$data[[region_key]])) {
      bundle_cache$data[[region_key]] <- load_region_bundle(
        data_dir = data_dir,
        region_id = input$region_id
      )
    }

    current_bundle(bundle_cache$data[[region_key]])
    selected_ids(character())
    interventions(empty_interventions())
    scenario_output(NULL)
  }, ignoreNULL = FALSE)

  observeEvent(current_bundle(), {
    bundle <- req(current_bundle())

    updateSelectizeInput(
      session = session,
      inputId = "tract_table_geoid",
      choices = set_names(bundle$tract_data$GEOID, bundle$tract_data$GEOID),
      selected = character(0),
      server = TRUE
    )
  })

  output$tract_table <- renderTable({
    bundle <- req(current_bundle())

    bundle$tract_data |>
      st_drop_geometry() |>
      select(GEOID, CountyName, res_obs, emp_obs, Q_obs, a_bl, b_bl, varphi_bl) |>
      slice_head(n = 15)
  })

  observeEvent(input$baseline_map_shape_click, {
    clicked_id <- input$baseline_map_shape_click$id
    current_ids <- selected_ids()

    if (clicked_id %in% current_ids) {
      selected_ids(setdiff(current_ids, clicked_id))
    } else {
      selected_ids(c(current_ids, clicked_id) |> unique() |> sort())
    }
  })

  observeEvent(input$add_table_selection, {
    table_ids <- input$tract_table_geoid

    if (length(table_ids) == 0) {
      return(invisible(NULL))
    }

    selected_ids(c(selected_ids(), table_ids) |> unique() |> sort())
  })

  observeEvent(input$clear_selection, {
    selected_ids(character())
  })

  observeEvent(input$add_intervention, {
    if (length(selected_ids()) == 0) {
      showNotification("Select at least one tract before adding an intervention.", type = "warning")
      return(invisible(NULL))
    }

    current_stack <- interventions()

    intervention_record <- build_intervention_record(
      tract_ids = selected_ids(),
      tract_set_name = input$tract_set_name,
      target_variable = input$target_variable,
      shock_percent = input$shock_percent,
      order = nrow(current_stack) + 1
    )

    interventions(bind_rows(current_stack, intervention_record))
    showNotification("Intervention added to the current scenario.", type = "message")
  })

  observeEvent(input$remove_last_intervention, {
    current_stack <- interventions()

    if (nrow(current_stack) == 0) {
      return(invisible(NULL))
    }

    interventions(slice_head(current_stack, n = nrow(current_stack) - 1))
  })

  observeEvent(input$clear_interventions, {
    interventions(empty_interventions())
  })

  output$selected_tract_summary <- renderPrint({
    ids <- selected_ids()

    if (length(ids) == 0) {
      cat("No tracts selected.")
    } else {
      cat(length(ids), "tract(s) selected\n")
      print(ids)
    }
  })

  output$intervention_table <- renderTable({
    interventions() |>
      mutate(
        shock_label = paste0(shock_percent, "%"),
        tract_ids = map_chr(tract_ids, ~ paste(.x, collapse = ", "))
      ) |>
      select(order, tract_set_name, target_variable, shock_label, tract_count)
  })

  output$intervention_count <- renderText({
    paste(nrow(interventions()), "intervention(s) added")
  })

  output$baseline_map <- renderLeaflet({
    bundle <- req(current_bundle())
    value_col <- req(input$baseline_map_var)

    map_data <- bundle$tract_data |>
      mutate(
        fill_value = .data[[value_col]],
        is_selected = GEOID %in% selected_ids(),
        popup_text = paste0(
          "<strong>GEOID:</strong> ", GEOID, "<br/>",
          "<strong>County:</strong> ", CountyName, "<br/>",
          "<strong>Value:</strong> ", round(fill_value, 3)
        )
      )

    pal <- build_palette(map_data$fill_value)

    leaflet(map_data) |>
      addProviderTiles("OpenStreetMap.Mapnik") |>
      addPolygons(
        layerId = ~GEOID,
        fillColor = ~pal(fill_value),
        fillOpacity = 0.8,
        color = ~if_else(is_selected, "#E07A00", "#4d4d4d"),
        weight = ~if_else(is_selected, 3, 1),
        popup = ~popup_text
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = map_data$fill_value,
        title = names(baseline_choices[baseline_choices == value_col])
      )
  })

  observeEvent(input$run_scenario, {
    bundle <- req(current_bundle())

    withProgress(message = "Running scenario", value = 0, {
      incProgress(0.2, detail = "Applying interventions")
      current_interventions <- interventions()

      incProgress(0.6, detail = "Solving counterfactual equilibrium")
      scenario_result <- run_counterfactual(bundle, current_interventions)

      incProgress(1, detail = "Preparing outputs")
      scenario_output(scenario_result)
      updateTabsetPanel(session, "main_view", selected = "results")
      updateTabsetPanel(session, "results_view", selected = "map")
    })
  })

  output$summary_table <- renderTable({
    validate(need(scenario_output(), "Run a scenario to view its summary."))
    scenario_output()$summary_table |>
      mutate(
        baseline = round(baseline, 4),
        counterfactual = round(counterfactual, 4),
        absolute_change = round(absolute_change, 4),
        percent_change = round(percent_change, 4)
      )
  })

  output$results_map <- renderLeaflet({
    validate(need(scenario_output(), "Run a scenario to view its results map."))

    result_sf <- scenario_output()$result_sf
    display_col <- resolve_result_column(
      base_name = input$result_map_var,
      display_mode = input$result_display_mode
    )

    map_data <- result_sf |>
      mutate(
        is_selected = GEOID %in% selected_ids(),
        fill_value = .data[[display_col]],
        fill_value = if_else(is.finite(fill_value), fill_value, NA_real_),
        fill_value = if_else(input$hide_selected_in_results & is_selected, NA_real_, fill_value),
        popup_text = paste0(
          "<strong>GEOID:</strong> ", GEOID, "<br/>",
          "<strong>County:</strong> ", CountyName, "<br/>",
          "<strong>Display value:</strong> ", format_result_values(fill_value, input$result_display_mode)
        )
      )

    pal <- build_palette(map_data$fill_value)

    leaflet(map_data) |>
      addProviderTiles("OpenStreetMap.Mapnik") |>
      addPolygons(
        layerId = ~GEOID,
        fillColor = ~pal(fill_value),
        fillOpacity = ~if_else(is_selected & input$hide_selected_in_results, 0.05, 0.8),
        color = ~if_else(is_selected & input$hide_selected_in_results, "#bdbdbd", "#4d4d4d"),
        weight = ~if_else(is_selected & input$hide_selected_in_results, 0.5, 1),
        popup = ~popup_text
      ) |>
      add_result_legend(
        pal = pal,
        values = map_data$fill_value,
        display_mode = input$result_display_mode,
        title = paste(
          names(result_choices[result_choices == input$result_map_var]),
          names(display_choices[display_choices == input$result_display_mode])
        )
      )
  })
}

shinyApp(ui = ui, server = server)
