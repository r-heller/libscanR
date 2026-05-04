# Export module
mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Export",
      shiny::downloadButton(ns("dl_spectra"), "Spectra (CSV)",
                            class = "btn-primary w-100 mb-2"),
      shiny::downloadButton(ns("dl_peaks"), "Peaks (CSV)",
                            class = "btn-primary w-100 mb-2"),
      shiny::downloadButton(ns("dl_calibration"), "Calibration (CSV)",
                            class = "btn-primary w-100 mb-2"),
      shiny::downloadButton(ns("dl_info"), "Sample info (CSV)",
                            class = "btn-primary w-100 mb-2")
    ),
    bslib::card(
      bslib::card_header("Export summary"),
      shiny::uiOutput(ns("summary"))
    )
  )
}

mod_export_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    output$summary <- shiny::renderUI({
      data <- rv$processed_data %||% rv$raw_data
      if (is.null(data)) {
        return(shiny::div(class = "libscanR-status",
                          "Nothing to export. Load data in the Import tab first."))
      }
      txt <- if (inherits(data, "libs_dataset")) {
        sprintf("Ready: %d spectra, %d channels.",
                data$n_spectra, data$n_channels)
      } else {
        sprintf("Ready: spectrum with %d channels, %d shots.",
                data$n_channels, data$n_shots)
      }
      shiny::div(class = "libscanR-status", txt)
    })

    output$dl_spectra <- shiny::downloadHandler(
      filename = function() "libscanR_spectra.csv",
      content = function(file) {
        data <- rv$processed_data %||% rv$raw_data
        shiny::req(data)
        if (inherits(data, "libs_spectrum")) {
          libscanR::ls_write_csv(data, file)
        } else {
          # Flatten: one column per sample
          mat <- data$intensity_matrix
          df <- data.frame(wavelength_nm = data$wavelength)
          for (i in seq_len(nrow(mat))) {
            df[[as.character(data$sample_info$sample_id[i])]] <- mat[i, ]
          }
          utils::write.csv(df, file, row.names = FALSE)
        }
      }
    )
    output$dl_peaks <- shiny::downloadHandler(
      filename = function() "libscanR_peaks.csv",
      content = function(file) {
        shiny::req(rv$peaks)
        utils::write.csv(rv$peaks, file, row.names = FALSE)
      }
    )
    output$dl_calibration <- shiny::downloadHandler(
      filename = function() "libscanR_calibration.csv",
      content = function(file) {
        cals <- rv$calibrations
        shiny::req(length(cals) > 0)
        rows <- lapply(cals, function(c) {
          data.frame(
            element = c$element,
            wavelength_nm = c$wavelength_nm,
            method = c$method,
            r_squared = c$r_squared %||% NA,
            lod = c$lod %||% NA,
            loq = c$loq %||% NA,
            unit = c$unit,
            n_standards = c$n_standards
          )
        })
        df <- do.call(rbind, rows)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_info <- shiny::downloadHandler(
      filename = function() "libscanR_sample_info.csv",
      content = function(file) {
        data <- rv$processed_data %||% rv$raw_data
        shiny::req(inherits(data, "libs_dataset"))
        utils::write.csv(data$sample_info, file, row.names = FALSE)
      }
    )
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
