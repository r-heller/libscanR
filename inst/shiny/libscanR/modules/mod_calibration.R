# Calibration module
mod_calibration_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Calibration",
      shiny::selectInput(ns("element"), "Element",
                         choices = c("Ca", "Na", "K", "Mg", "Fe", "Zn",
                                     "Cu", "Pb", "Cd", "Cr"),
                         selected = "Ca"),
      shiny::numericInput(ns("line_nm"), "Emission line (nm)",
                          value = 393.37, step = 0.01),
      shiny::selectInput(ns("method"), "Method",
                         choices = c("univariate", "internal_std", "pls"),
                         selected = "univariate"),
      shiny::numericInput(ns("internal_std"), "Internal std line (nm)",
                          value = 589.00, step = 0.01),
      shiny::numericInput(ns("window"), "Integration window (nm)",
                          value = 1, step = 0.1),
      shiny::textInput(ns("conc_col"), "Concentration column in sample_info",
                       value = "concentration"),
      shiny::actionButton(ns("build"), "Build calibration",
                          class = "btn-primary w-100")
    ),
    bslib::card(
      bslib::card_header("Calibration curve"),
      shiny::plotOutput(ns("curve"), height = "400px")
    ),
    bslib::card(
      bslib::card_header("Figures of merit"),
      shiny::verbatimTextOutput(ns("stats"))
    )
  )
}

mod_calibration_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$build, {
      ds <- rv$processed_data %||% rv$raw_data
      shiny::req(ds, inherits(ds, "libs_dataset"))
      if (!(input$conc_col %in% names(ds$sample_info))) {
        shiny::showNotification(paste("Column '", input$conc_col,
                                      "' not found."), type = "error")
        return()
      }
      concs <- ds$sample_info[[input$conc_col]]
      cal <- tryCatch(
        libscanR::ls_calibrate(ds, input$element, input$line_nm,
                               concentrations = concs,
                               method = input$method,
                               internal_std_nm = input$internal_std,
                               window_nm = input$window,
                               verbose = FALSE),
        error = function(e) {
          shiny::showNotification(conditionMessage(e), type = "error")
          NULL
        })
      if (!is.null(cal)) {
        rv$calibrations[[input$element]] <- cal
      }
    })

    output$curve <- shiny::renderPlot({
      cal <- rv$calibrations[[input$element]]
      shiny::req(cal)
      libscanR::ls_plot_calibration(cal)
    })

    output$stats <- shiny::renderPrint({
      cal <- rv$calibrations[[input$element]]
      shiny::req(cal)
      summary(cal)
    })
  })
}
