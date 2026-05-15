app_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
app_dir <- if (!is.null(app_file)) {
  dirname(normalizePath(app_file, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
common_candidates <- unique(c(
  file.path(dirname(app_dir), "common.R"),
  file.path(getwd(), "interactive_graphics", "common.R")
))
common_path <- common_candidates[file.exists(common_candidates)][1]
if (is.na(common_path) || !nzchar(common_path)) {
  stop("Не удалось найти interactive_graphics/common.R.")
}
source(common_path, local = TRUE)

ui <- shiny::fluidPage(
  shiny::titlePanel("Зависимость ошибки от частоты"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::numericInput("freq_sigma", "Уровень шума", value = 0.2, min = 0, step = 0.01),
      shiny::numericInput("freq_n_rep", "Размер выборки", value = 1, min = 1, step = 1),
      shiny::numericInput(
        "freq_n_workers",
        "Число потоков",
        value = available_worker_count(),
        min = 1,
        step = 1
      ),
      shiny::numericInput("freq_n", "Размер вектора", value = 100, min = 8, step = 1),
      shiny::numericInput("freq_seed", "Seed", value = 1111, min = 1, step = 1),
      shiny::actionButton("run_frequency", "Построить графики")
    ),
    shiny::mainPanel(
      shiny::h4("Краткая сводка"),
      shiny::tableOutput("freq_head_table"),
      shiny::h4("Ошибка локализации от частоты"),
      shiny::plotOutput("freq_error_plot", height = "380px"),
      shiny::h4("Дополнительные метрики"),
      shiny::plotOutput("freq_metrics_plot", height = "760px"),
      shiny::h4("Частоты с наибольшей ошибкой"),
      shiny::tableOutput("freq_top_error_table")
    )
  )
)

server <- function(input, output, session) {
  frequency_result <- shiny::eventReactive(input$run_frequency, {
    shiny::withProgress(message = "Расчет частотного профиля", value = 0, {
      incProgress(0.2, detail = "Генерация выборки")
      samples <- run_unit_frequency_experiment(
        n = as.integer(input$freq_n),
        n_rep = as.integer(input$freq_n_rep),
        sigma = as.numeric(input$freq_sigma),
        seed = as.integer(input$freq_seed),
        n_workers = input$freq_n_workers
      )
      incProgress(0.8, detail = "Агрегация")
      list(
        samples = samples,
        summary = summarise_unit_frequency(samples)
      )
    })
  }, ignoreNULL = FALSE)

  output$freq_head_table <- shiny::renderTable({
    shiny::req(frequency_result())
    head(frequency_result()$summary, 10)
  }, digits = 6)

  output$freq_error_plot <- shiny::renderPlot({
    shiny::req(frequency_result())
    summary_df <- frequency_result()$summary
    plot(
      summary_df$omega,
      summary_df$argmax_error,
      type = "b",
      pch = 19,
      col = "#1f4e79",
      xlab = expression(omega[j] == 2*pi*(j-1)/N),
      ylab = expression(abs(hat(j) - j)),
      main = "Ошибка локализации максимума"
    )
    grid(col = "gray85")
    abline(stats::lm(argmax_error ~ omega, data = summary_df), col = "red", lwd = 2)
  })

  output$freq_metrics_plot <- shiny::renderPlot({
    shiny::req(frequency_result())
    summary_df <- frequency_result()$summary
    old_par <- par(mfrow = c(3, 1), mar = c(4, 4, 2.2, 1))
    on.exit(par(old_par), add = TRUE)

    plot(
      summary_df$omega,
      summary_df$mse,
      type = "b",
      pch = 19,
      col = "#2a6f37",
      xlab = expression(omega[j]),
      ylab = "MSE",
      main = "MSE восстановления"
    )
    grid(col = "gray85")

    plot(
      summary_df$omega,
      summary_df$near_mass_share,
      type = "b",
      pch = 19,
      col = "#7b3f00",
      xlab = expression(omega[j]),
      ylab = "Доля массы",
      main = "Доля положительной массы в окрестности j±1"
    )
    grid(col = "gray85")

    plot(
      summary_df$omega,
      summary_df$active_share_0001,
      type = "b",
      pch = 19,
      col = "#6b3fa0",
      xlab = expression(omega[j]),
      ylab = "Доля активных",
      main = expression(hat(x) > 10^{-3})
    )
    grid(col = "gray85")
  })

  output$freq_top_error_table <- shiny::renderTable({
    shiny::req(frequency_result())
    summary_df <- frequency_result()$summary
    summary_df[order(-summary_df$argmax_error), ][1:min(10, nrow(summary_df)), ]
  }, digits = 6)
}

shiny::shinyApp(ui = ui, server = server)
