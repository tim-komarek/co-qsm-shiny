library(shiny)
library(leaflet)
library(tidyverse)
library(sf)

app_dir <- normalizePath(getwd())
data_dir <- file.path(app_dir, "data")

source(file.path(app_dir, "R", "model_bundle_helpers.R"))
source(file.path(app_dir, "R", "scenario_helpers.R"))
source_model_functions()

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
  leaflet::colorNumeric(
    palette = "viridis",
    domain = values,
    na.color = "#d9d9d9"
  )
}

ui <- fluidPage(
  titlePanel("Colorado Quantitative Spatial Model"),
  if (nrow(bundle_index) == 0) {
    fluidRow(
      column(
        width = 12,
        h4("No Shiny region bundles found."),
        p("Run `Rscript Shiny/build_region_bundles.R` from the project root, then restart the app.")
      )
    )
  } else {
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = "region_id",
          label = "Region",
          choices = set_names(bundle_index$region_id, bundle_index$region_name),
          selected = bundle_index$region_id[[1]]
        ),
        tags$hr(),
        h4("Baseline Map"),
        selectInput(
          inputId = "baseline_map_var",
          label = "Baseline variable",
          choices = baseline_choices,
          selected = "res_obs"
        ),
        actionButton("clear_selection", "Clear selected tracts"),
        tags$hr(),
        h4("Intervention Builder"),
        textInput(
          inputId = "tract_set_name",
          label = "Tract set name",
          value = "selected_tracts"
        ),
        selectInput(
          inputId = "target_variable",
          label = "Policy variable",
          choices = policy_choices,
          selected = "b"
        ),
        sliderInput(
          inputId = "shock_percent",
          label = "Percent change",
          min = -95,
          max = 100,
          value = 10,
          step = 1,
          post = "%"
        ),
        actionButton("add_table_selection", "Add table selection to current tract set"),
        actionButton("add_intervention", "Add intervention"),
        actionButton("remove_last_intervention", "Remove last intervention"),
        actionButton("clear_interventions", "Clear interventions"),
        tags$hr(),
        actionButton("run_scenario", "Run scenario", class = "btn-primary")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Baseline Explorer",
            br(),
            fluidRow(
              column(
                width = 8,
                leafletOutput("baseline_map", height = 700)
              ),
              column(
                width = 4,
                h4("Selected Tracts"),
                verbatimTextOutput("selected_tract_summary"),
                selectizeInput(
                  inputId = "tract_table_geoid",
                  label = "Add tracts from list",
                  choices = NULL,
                  multiple = TRUE
                ),
                h4("Tract Table"),
                tableOutput("tract_table")
              )
            )
          ),
          tabPanel(
            "Scenario Results",
            br(),
            fluidRow(
              column(
                width = 3,
                selectInput(
                  inputId = "result_map_var",
                  label = "Result variable",
                  choices = result_choices,
                  selected = "u"
                ),
                selectInput(
                  inputId = "result_display_mode",
                  label = "Display mode",
                  choices = display_choices,
                  selected = "pct"
                ),
                checkboxInput(
                  inputId = "hide_selected_in_results",
                  label = "Hide selected tracts in result map",
                  value = FALSE
                ),
                h4("Intervention Stack"),
                tableOutput("intervention_table"),
                h4("Scenario Summary"),
                tableOutput("summary_table")
              ),
              column(
                width = 9,
                leafletOutput("results_map", height = 700)
              )
            )
          )
        )
      )
    )
  }
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
      addProviderTiles("CartoDB.Positron") |>
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
    })
  })

  output$summary_table <- renderTable({
    req(scenario_output())
    scenario_output()$summary_table |>
      mutate(
        baseline = round(baseline, 4),
        counterfactual = round(counterfactual, 4),
        absolute_change = round(absolute_change, 4),
        percent_change = round(percent_change, 4)
      )
  })

  output$results_map <- renderLeaflet({
    req(scenario_output())

    result_sf <- scenario_output()$result_sf
    display_col <- resolve_result_column(
      base_name = input$result_map_var,
      display_mode = input$result_display_mode
    )

    map_data <- result_sf |>
      mutate(
        is_selected = GEOID %in% selected_ids(),
        fill_value = .data[[display_col]],
        fill_value = if_else(input$hide_selected_in_results & is_selected, NA_real_, fill_value),
        popup_text = paste0(
          "<strong>GEOID:</strong> ", GEOID, "<br/>",
          "<strong>County:</strong> ", CountyName, "<br/>",
          "<strong>Display value:</strong> ", round(fill_value, 3)
        )
      )

    pal <- build_palette(map_data$fill_value)

    leaflet(map_data) |>
      addProviderTiles("CartoDB.Positron") |>
      addPolygons(
        layerId = ~GEOID,
        fillColor = ~pal(fill_value),
        fillOpacity = ~if_else(is_selected & input$hide_selected_in_results, 0.05, 0.8),
        color = ~if_else(is_selected & input$hide_selected_in_results, "#bdbdbd", "#4d4d4d"),
        weight = ~if_else(is_selected & input$hide_selected_in_results, 0.5, 1),
        popup = ~popup_text
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = map_data$fill_value,
        title = paste(
          names(result_choices[result_choices == input$result_map_var]),
          names(display_choices[display_choices == input$result_display_mode])
        )
      )
  })
}

shinyApp(ui = ui, server = server)
