# Peaks module
mod_peaks_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Peak detection",
      shiny::sliderInput(ns("snr"), "SNR threshold",
                         min = 1, max = 20, value = 5, step = 0.5),
      shiny::sliderInput(ns("prom"), "Min prominence (fraction of max)",
                         min = 0.001, max = 0.2, value = 0.01, step = 0.001),
      shiny::sliderInput(ns("dist"), "Min distance (nm)",
                         min = 0.1, max = 3, value = 0.5, step = 0.1),
      shiny::selectInput(ns("elements"), "Annotate elements (NIST)",
                         choices = c("Ca", "Na", "K", "Mg", "Fe",
                                     "Zn", "Cu", "P", "C", "H", "O",
                                     "N", "Mn", "Sr", "Ba", "Pb", "Cd",
                                     "Cr", "Al", "Si", "Li", "Ti", "Cl"),
                         multiple = TRUE,
                         selected = c("Ca", "Na", "K", "Fe")),
      shiny::actionButton(ns("detect"), "Detect peaks",
                          class = "btn-primary w-100")
    ),
    bslib::card(
      bslib::card_header("Spectrum with peaks"),
      plotly::plotlyOutput(ns("plot"), height = "400px")
    ),
    bslib::card(
      bslib::card_header("Identified peaks"),
      DT::DTOutput(ns("table"))
    )
  )
}

mod_peaks_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$detect, {
      spec <- rv$processed_data %||% rv$raw_data
      shiny::req(spec)
      if (inherits(spec, "libs_dataset")) spec <- spec$spectra[[1]]
      pk <- libscanR::ls_find_peaks(spec,
                                    snr_threshold = input$snr,
                                    min_prominence = input$prom,
                                    min_distance_nm = input$dist)
      if (nrow(pk) > 0 && length(input$elements) > 0) {
        pk <- libscanR::ls_identify_peaks(pk, elements = input$elements)
      }
      rv$peaks <- pk
    })

    output$plot <- plotly::renderPlotly({
      spec <- rv$processed_data %||% rv$raw_data
      shiny::req(spec)
      if (inherits(spec, "libs_dataset")) spec <- spec$spectra[[1]]
      p <- libscanR::ls_plot_spectrum(
        spec,
        show_peaks = !is.null(rv$peaks) && nrow(rv$peaks) > 0,
        snr_threshold = input$snr,
        show_elements = input$elements
      )
      plotly::ggplotly(p)
    })

    output$table <- DT::renderDT({
      shiny::req(rv$peaks)
      df <- rv$peaks
      df[] <- lapply(df, function(col) {
        if (is.numeric(col)) round(col, 4) else col
      })
      DT::datatable(df, options = list(pageLength = 15, dom = "tip"))
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
