# Chemometrics module
mod_chemometrics_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_tab(
    id = ns("tabs"),
    bslib::nav_panel("PCA",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "PCA",
          shiny::numericInput(ns("pca_n"), "Components", value = 5,
                              min = 2, max = 20),
          shiny::textInput(ns("pca_color"), "Color by (column)",
                           value = "tissue"),
          shiny::actionButton(ns("run_pca"), "Run PCA",
                              class = "btn-primary w-100")
        ),
        bslib::card(bslib::card_header("Scores"),
                    plotly::plotlyOutput(ns("pca_scores"),
                                         height = "400px")),
        bslib::card(bslib::card_header("Scree"),
                    shiny::plotOutput(ns("pca_scree"), height = "250px")),
        bslib::card(bslib::card_header("Loadings (PC1)"),
                    shiny::plotOutput(ns("pca_loadings"), height = "300px"))
      )
    ),
    bslib::nav_panel("PLS-DA",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "PLS-DA",
          shiny::textInput(ns("plsda_col"), "Class column",
                           value = "tissue"),
          shiny::numericInput(ns("plsda_n"), "Components",
                              value = 3, min = 2, max = 20),
          shiny::actionButton(ns("run_plsda"), "Train",
                              class = "btn-primary w-100")
        ),
        bslib::card(bslib::card_header("Scores"),
                    shiny::plotOutput(ns("plsda_scores"), height = "400px")),
        bslib::card(bslib::card_header("Confusion matrix"),
                    shiny::plotOutput(ns("plsda_cm"), height = "300px"))
      )
    ),
    bslib::nav_panel("Cluster",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Clustering",
          shiny::selectInput(ns("cluster_method"), "Method",
                             choices = c("kmeans", "hclust")),
          shiny::sliderInput(ns("cluster_k"), "k", min = 2, max = 15,
                             value = 3),
          shiny::actionButton(ns("run_cluster"), "Run",
                              class = "btn-primary w-100")
        ),
        bslib::card(bslib::card_header("Cluster assignments"),
                    DT::DTOutput(ns("cluster_tbl")))
      )
    )
  )
}

mod_chemometrics_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$run_pca, {
      ds <- rv$processed_data %||% rv$raw_data
      shiny::req(ds, inherits(ds, "libs_dataset"))
      rv$pca_result <- libscanR::ls_pca(ds, n_components = input$pca_n)
    })

    output$pca_scores <- plotly::renderPlotly({
      shiny::req(rv$pca_result)
      cb <- input$pca_color
      if (!cb %in% names(rv$pca_result$sample_info)) cb <- NULL
      p <- libscanR::ls_plot_pca(rv$pca_result, color_by = cb)
      plotly::ggplotly(p)
    })
    output$pca_scree <- shiny::renderPlot({
      shiny::req(rv$pca_result)
      libscanR::ls_plot_scree(rv$pca_result)
    })
    output$pca_loadings <- shiny::renderPlot({
      shiny::req(rv$pca_result)
      libscanR::ls_plot_loadings(rv$pca_result, pc = 1, n_top = 5)
    })

    shiny::observeEvent(input$run_plsda, {
      ds <- rv$processed_data %||% rv$raw_data
      shiny::req(ds, inherits(ds, "libs_dataset"))
      plsda <- tryCatch(
        libscanR::ls_plsda(ds, grouping = input$plsda_col,
                           n_components = input$plsda_n,
                           validation = "none"),
        error = function(e) {
          shiny::showNotification(conditionMessage(e), type = "error")
          NULL
        })
      rv$plsda_result <- plsda
    })

    output$plsda_scores <- shiny::renderPlot({
      shiny::req(rv$plsda_result)
      libscanR::ls_plot_plsda(rv$plsda_result, type = "scores")
    })
    output$plsda_cm <- shiny::renderPlot({
      shiny::req(rv$plsda_result)
      libscanR::ls_plot_plsda(rv$plsda_result, type = "confusion")
    })

    shiny::observeEvent(input$run_cluster, {
      ds <- rv$processed_data %||% rv$raw_data
      shiny::req(ds, inherits(ds, "libs_dataset"))
      rv$cluster_result <- libscanR::ls_cluster(ds,
                                                method = input$cluster_method,
                                                k = input$cluster_k)
    })

    output$cluster_tbl <- DT::renderDT({
      shiny::req(rv$cluster_result)
      df <- data.frame(
        sample_id = rv$cluster_result$sample_info$sample_id,
        cluster = rv$cluster_result$cluster
      )
      DT::datatable(df, options = list(pageLength = 15, dom = "tip"))
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
