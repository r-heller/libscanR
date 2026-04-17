# Import module
mod_import_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Data import",
      shiny::selectInput(ns("format"), "Format",
                         choices = c("Auto-detect" = "auto",
                                     "Generic CSV/TSV" = "generic",
                                     "SciAps" = "sciaps",
                                     "Applied Spectra" = "aurora")),
      shiny::fileInput(ns("file"), "Upload spectrum (CSV/TSV/TXT)",
                       accept = c(".csv", ".tsv", ".txt"),
                       multiple = TRUE),
      shiny::hr(),
      shiny::strong("Or load example data:"),
      shiny::actionButton(ns("load_tissue"), "Tissue dataset",
                          class = "btn-outline-primary w-100"),
      shiny::actionButton(ns("load_cal"), "Calibration dataset",
                          class = "btn-outline-primary w-100 mt-2"),
      shiny::actionButton(ns("load_spatial"), "Spatial dataset",
                          class = "btn-outline-primary w-100 mt-2")
    ),
    bslib::card(
      bslib::card_header("Import status"),
      shiny::uiOutput(ns("status"))
    ),
    bslib::card(
      bslib::card_header("Preview"),
      plotly::plotlyOutput(ns("preview"), height = "400px")
    ),
    bslib::card(
      bslib::card_header("Sample info"),
      DT::DTOutput(ns("info"))
    )
  )
}

mod_import_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$file, {
      shiny::req(input$file)
      specs <- list()
      for (i in seq_len(nrow(input$file))) {
        path <- input$file$datapath[i]
        spec <- tryCatch({
          switch(input$format,
                 auto = libscanR::ls_read_auto(path, verbose = FALSE),
                 generic = libscanR::ls_read_spectrum(path, verbose = FALSE),
                 sciaps = libscanR::ls_read_sciaps(path, verbose = FALSE),
                 aurora = libscanR::ls_read_aurora(path, verbose = FALSE))
        }, error = function(e) NULL)
        if (!is.null(spec)) specs[[length(specs) + 1]] <- spec
      }
      if (length(specs) == 0) return()
      if (length(specs) == 1) {
        rv$raw_data <- specs[[1]]
      } else {
        rv$raw_data <- libscanR::ls_dataset(specs)
      }
    })

    shiny::observeEvent(input$load_tissue, {
      rv$raw_data <- libscanR::ls_example_data("tissue")
    })
    shiny::observeEvent(input$load_cal, {
      rv$raw_data <- libscanR::ls_example_data("calibration")
    })
    shiny::observeEvent(input$load_spatial, {
      rv$raw_data <- libscanR::ls_example_data("spatial")
    })

    output$status <- shiny::renderUI({
      if (is.null(rv$raw_data)) {
        return(shiny::div(class = "libscanR-status",
                          "No data loaded yet. Upload a spectrum or load an example."))
      }
      if (inherits(rv$raw_data, "libs_spectrum")) {
        txt <- sprintf("Spectrum loaded: %d channels, %d shots, range %.1f-%.1f nm.",
                       rv$raw_data$n_channels, rv$raw_data$n_shots,
                       rv$raw_data$range_nm[1], rv$raw_data$range_nm[2])
      } else {
        txt <- sprintf("Dataset loaded: %d spectra, %d channels, range %.1f-%.1f nm.",
                       rv$raw_data$n_spectra, rv$raw_data$n_channels,
                       rv$raw_data$range_nm[1], rv$raw_data$range_nm[2])
      }
      shiny::div(class = "libscanR-status", txt)
    })

    output$preview <- plotly::renderPlotly({
      shiny::req(rv$raw_data)
      if (inherits(rv$raw_data, "libs_spectrum")) {
        p <- libscanR::ls_plot_spectrum(rv$raw_data)
      } else {
        n <- min(10, rv$raw_data$n_spectra)
        p <- libscanR::ls_plot_overlay(rv$raw_data[seq_len(n)],
                                       color_by = NULL)
      }
      plotly::ggplotly(p)
    })

    output$info <- DT::renderDT({
      shiny::req(rv$raw_data)
      if (inherits(rv$raw_data, "libs_spectrum")) {
        df <- data.frame(
          key = names(rv$raw_data$metadata),
          value = vapply(rv$raw_data$metadata, function(v) {
            paste(as.character(v), collapse = ", ")
          }, character(1)),
          stringsAsFactors = FALSE
        )
      } else {
        df <- as.data.frame(rv$raw_data$sample_info)
      }
      DT::datatable(df, options = list(pageLength = 10, dom = "tip"))
    })
  })
}
