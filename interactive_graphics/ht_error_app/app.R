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

method_choices <- detector_labels(available_detectors(num_of_lines = 1L, include_esprit = FALSE))

resolve_cssa_method_names <- function(method_names, use_esprit = FALSE) {
  method_names <- as.character(method_names)
  if (!isTRUE(use_esprit)) {
    return(method_names)
  }

  method_names[method_names == "cssa_row_row"] <- "cssa_row_row_esprit"
  method_names[method_names == "cssa_col_row"] <- "cssa_col_row_esprit"
  method_names
}

ui <- shiny::fluidPage(
  shiny::titlePanel("HT: сравнение методов и поиск больших ошибок"),
  shiny::tabsetPanel(
    shiny::tabPanel(
      "Сравнение методов",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::numericInput("exp_n_rep", "Размер выборки", value = 5, min = 1, step = 1),
          shiny::numericInput(
            "exp_n_workers",
            "Число потоков",
            value = available_worker_count(),
            min = 1,
            step = 1
          ),
          shiny::numericInput("exp_sigma", "Уровень шума", value = 0.2, min = 0, step = 0.01),
          shiny::numericInput("exp_n_row", "Число строк матрицы", value = 100, min = 10, step = 1),
          shiny::numericInput("exp_n_col", "Число столбцов матрицы", value = 100, min = 10, step = 1),
          shiny::numericInput("exp_num_lines", "Количество прямых", value = 1, min = 1, step = 1),
          shiny::checkboxInput("exp_cssa_esprit", "CSSA: ESPRIT + corrected frequency", value = FALSE),
          shiny::numericInput("exp_max_row_k", "MAX.row: points per row", value = 1, min = 1, step = 1),
          shiny::checkboxInput(
            "exp_use_standard_single",
            "Для одной прямой использовать 4 стандартные конфигурации",
            value = TRUE
          ),
          shiny::conditionalPanel(
            "input.exp_num_lines > 1 || !input.exp_use_standard_single",
            shiny::textAreaInput(
              "exp_line_text",
              "Параметры прямых: a,b[,intensity], по одной на строку",
              value = "1,0\n-1,101\n2,-1\n2,-40",
              rows = 6
            )
          ),
          shiny::selectInput(
            "exp_ht_type",
            "Метод HT",
            choices = c("Обычный" = "ordinary", "Weighted" = "weighted"),
            selected = "ordinary"
          ),
          shiny::numericInput("exp_rho_step", "Шаг rho", value = 1, min = 0.001, step = 0.1),
          shiny::numericInput("exp_theta_step", "Шаг theta", value = 0.01, min = 0.0001, step = 0.001),
          shiny::selectInput(
            "exp_threshold_mode",
            "Порог",
            choices = c(
              "Авто: SD остатка row.row" = "auto_rowrow_sd",
              "Без порога" = "none",
              "Вручную" = "manual",
              "p-квантиль" = "quantile"
            ),
            selected = "none"
          ),
          shiny::conditionalPanel(
            "input.exp_threshold_mode == 'manual'",
            shiny::numericInput("exp_threshold_value", "Ручной порог", value = 0.1, min = 0, step = 0.01)
          ),
          shiny::conditionalPanel(
            "input.exp_threshold_mode == 'quantile'",
            shiny::numericInput("exp_threshold_quantile_p", "p для квантиля", value = 0.9, min = 0, max = 1, step = 0.01)
          ),
          shiny::conditionalPanel(
            "input.exp_threshold_mode == 'auto_rowrow_sd'",
            shiny::numericInput(
              "exp_threshold_multiplier",
              "Множитель для авто-порога",
              value = 1,
              min = 0,
              step = 0.1
            )
          ),
          shiny::conditionalPanel(
            "input.exp_threshold_mode != 'none'",
            shiny::checkboxInput("exp_threshold_all_methods", "Выбрать все", value = TRUE),
            shiny::checkboxGroupInput(
              "exp_threshold_methods",
              "Применить порог к методам",
              choices = method_choices,
              selected = unname(method_choices)
            )
          ),
          shiny::checkboxGroupInput(
            "exp_methods",
            "Методы обработки",
            choices = method_choices,
            selected = unname(method_choices)
          ),
          shiny::helpText(if (!bm3d_available) "BM3D недоступен в текущем Python-окружении." else NULL),
          shiny::actionButton("run_experiment", "Запустить эксперимент")
        ),
        shiny::mainPanel(
          shiny::h4("Предпросмотр конфигураций"),
          shiny::verbatimTextOutput("exp_config_message"),
          shiny::plotOutput("exp_config_preview", height = "420px"),
          shiny::hr(),
          shiny::h4("Сводка по ошибкам"),
          shiny::tableOutput("exp_summary_table"),
          shiny::hr(),
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Mean dr",
              shiny::tags$p(shiny::strong("Таблица ошибок")),
              shiny::tableOutput("exp_mean_dr"),
              shiny::tags$p(shiny::strong("Таблица ошибок, деленных на минимумы по строкам")),
              shiny::tableOutput("exp_mean_dr_norm")
            ),
            shiny::tabPanel(
              "Median dr",
              shiny::tags$p(shiny::strong("Таблица ошибок")),
              shiny::tableOutput("exp_median_dr"),
              shiny::tags$p(shiny::strong("Таблица ошибок, деленных на минимумы по строкам")),
              shiny::tableOutput("exp_median_dr_norm")
            ),
            shiny::tabPanel(
              "Mean dtheta",
              shiny::tags$p(shiny::strong("Таблица ошибок")),
              shiny::tableOutput("exp_mean_dtheta"),
              shiny::tags$p(shiny::strong("Таблица ошибок, деленных на минимумы по строкам")),
              shiny::tableOutput("exp_mean_dtheta_norm")
            ),
            shiny::tabPanel(
              "Median dtheta",
              shiny::tags$p(shiny::strong("Таблица ошибок")),
              shiny::tableOutput("exp_median_dtheta"),
              shiny::tags$p(shiny::strong("Таблица ошибок, деленных на минимумы по строкам")),
              shiny::tableOutput("exp_median_dtheta_norm")
            ),
            shiny::tabPanel("Активные пиксели", shiny::tableOutput("exp_active_pixels"), shiny::tableOutput("exp_active_weight"))
          )
        )
      )
    ),
    shiny::tabPanel(
      "Поиск больших ошибок",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::numericInput("find_sigma", "Уровень шума", value = 0.2, min = 0, step = 0.01),
          shiny::numericInput("find_n_row", "Число строк матрицы", value = 100, min = 10, step = 1),
          shiny::numericInput("find_n_col", "Число столбцов матрицы", value = 100, min = 10, step = 1),
          shiny::selectInput(
            "find_config_source",
            "Конфигурация",
            choices = c("Стандартная (1 прямая)" = "standard", "Пользовательская" = "custom"),
            selected = "standard"
          ),
          shiny::conditionalPanel(
            "input.find_config_source == 'standard'",
            shiny::selectInput(
              "find_standard_config",
              "Стандартная прямая",
              choices = c(
                "diag_pos_full" = "diag_pos_full",
                "diag_neg_full" = "diag_neg_full",
                "steep_full" = "steep_full",
                "steep_shifted" = "steep_shifted"
              )
            )
          ),
          shiny::conditionalPanel(
            "input.find_config_source == 'custom'",
            shiny::numericInput("find_num_lines", "Количество прямых", value = 1, min = 1, step = 1),
            shiny::textAreaInput(
              "find_line_text",
              "Параметры прямых: a,b[,intensity]",
              value = "2,-1",
              rows = 5
            )
          ),
          shiny::selectInput("find_method", "Метод обработки", choices = method_choices),
          shiny::checkboxInput("find_cssa_esprit", "CSSA: ESPRIT + corrected frequency", value = FALSE),
          shiny::numericInput("find_max_row_k", "MAX.row: points per row", value = 1, min = 1, step = 1),
          shiny::selectInput(
            "find_ht_type",
            "Метод HT",
            choices = c("Обычный" = "ordinary", "Weighted" = "weighted"),
            selected = "ordinary"
          ),
          shiny::numericInput("find_rho_step", "Шаг rho", value = 1, min = 0.001, step = 0.1),
          shiny::numericInput("find_theta_step", "Шаг theta", value = 0.01, min = 0.0001, step = 0.001),
          shiny::selectInput(
            "find_threshold_mode",
            "Порог",
            choices = c(
              "Авто: SD остатка row.row" = "auto_rowrow_sd",
              "Без порога" = "none",
              "Вручную" = "manual",
              "p-квантиль" = "quantile"
            ),
            selected = "auto_rowrow_sd"
          ),
          shiny::conditionalPanel(
            "input.find_threshold_mode == 'manual'",
            shiny::numericInput("find_threshold_value", "Ручной порог", value = 0.1, min = 0, step = 0.01)
          ),
          shiny::conditionalPanel(
            "input.find_threshold_mode == 'quantile'",
            shiny::numericInput("find_threshold_quantile_p", "p для квантиля", value = 0.9, min = 0, max = 1, step = 0.01)
          ),
          shiny::conditionalPanel(
            "input.find_threshold_mode == 'auto_rowrow_sd'",
            shiny::numericInput(
              "find_threshold_multiplier",
              "Множитель для авто-порога",
              value = 1,
              min = 0,
              step = 0.1
            )
          ),
          shiny::numericInput("find_num_maxima", "Количество максимумов", value = 1, min = 1, step = 1),
          shiny::numericInput("find_n_search", "Число шумовых реализаций", value = 100, min = 1, step = 1),
          shiny::numericInput(
            "find_n_workers",
            "Число потоков",
            value = available_worker_count(),
            min = 1,
            step = 1
          ),
          shiny::numericInput("find_factor", "Порог по отношению к идеалу", value = 5, min = 1, step = 0.5),
          shiny::numericInput("find_max_cases", "Максимум найденных случаев", value = 5, min = 1, step = 1),
          shiny::actionButton("run_find", "Искать большие ошибки")
        ),
        shiny::mainPanel(
          shiny::h4("Идеальная ошибка из-за дискретизации"),
          shiny::tableOutput("find_ideal_table"),
          shiny::hr(),
          shiny::h4("Найденные случаи"),
          shiny::tableOutput("find_cases_table"),
          shiny::uiOutput("find_case_selector"),
          shiny::tableOutput("find_selected_case_meta"),
          shiny::plotOutput("find_processed_plot", height = "420px"),
          shiny::plotOutput("find_accumulator_plot", height = "420px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  shiny::observeEvent(input$exp_methods, {
    selected_methods <- input$exp_methods %||% character(0)
    cur_threshold_methods <- input$exp_threshold_methods %||% character(0)
    keep_selected <- intersect(cur_threshold_methods, selected_methods)
    if (isTRUE(input$exp_threshold_all_methods)) {
      keep_selected <- selected_methods
    }
    shiny::updateCheckboxGroupInput(
      session,
      "exp_threshold_methods",
      choices = method_choices[unname(method_choices) %in% selected_methods],
      selected = keep_selected
    )
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$exp_threshold_all_methods, {
    selected_methods <- input$exp_methods %||% character(0)
    shiny::updateCheckboxGroupInput(
      session,
      "exp_threshold_methods",
      selected = if (isTRUE(input$exp_threshold_all_methods)) selected_methods else character(0)
    )
  }, ignoreInit = TRUE)

  experiment_config_list <- shiny::reactive({
    build_config_list(
      num_lines = input$exp_num_lines,
      n_row = input$exp_n_row,
      n_col = input$exp_n_col,
      use_standard_single = input$exp_use_standard_single,
      line_text = input$exp_line_text
    )
  })

  output$exp_config_message <- shiny::renderText({
    cfgs <- tryCatch(experiment_config_list(), error = function(e) e)
    if (inherits(cfgs, "error")) {
      return(paste("Ошибка конфигурации:", cfgs$message))
    }
    paste("Конфигураций:", length(cfgs))
  })

  output$exp_config_preview <- shiny::renderPlot({
    cfgs <- tryCatch(experiment_config_list(), error = function(e) e)
    if (inherits(cfgs, "error")) {
      plot.new()
      text(0.5, 0.5, cfgs$message)
      return(invisible(NULL))
    }
    draw_config_preview(cfgs)
  })

  experiment_result <- shiny::eventReactive(input$run_experiment, {
    cfgs <- experiment_config_list()
    chosen_labels <- input$exp_methods
    method_map <- setNames(names(method_choices), unname(method_choices))
    method_names <- unname(method_map[chosen_labels])
    threshold_method_names <- unname(method_map[input$exp_threshold_methods %||% character(0)])
    threshold_method_names <- threshold_method_names[!is.na(threshold_method_names)]
    method_names <- resolve_cssa_method_names(method_names, input$exp_cssa_esprit)
    threshold_method_names <- resolve_cssa_method_names(threshold_method_names, input$exp_cssa_esprit)

    shiny::withProgress(message = "Запуск эксперимента", value = 0, {
      incProgress(0.1, detail = "Подготовка")
      samples <- run_ht_experiment_app(
        config_list = cfgs,
        method_names = method_names,
        sigma_noise = input$exp_sigma,
        num_of_lines = input$exp_num_lines,
        rho_step_ht = input$exp_rho_step,
        theta_step_ht = input$exp_theta_step,
        n_rep = input$exp_n_rep,
        ht_type = input$exp_ht_type,
        threshold_mode = input$exp_threshold_mode,
        threshold_value = input$exp_threshold_value,
        threshold_quantile_p = input$exp_threshold_quantile_p,
        threshold_multiplier = input$exp_threshold_multiplier,
        threshold_method_names = threshold_method_names,
        n_workers = input$exp_n_workers,
        max_row_k = input$exp_max_row_k
      )
      incProgress(0.9, detail = "Сводка")
      list(
        samples = samples,
        summary = summarise_ht_errors_app(samples)
      )
    })
  })

  output$exp_summary_table <- shiny::renderTable({
    shiny::req(experiment_result())
    experiment_result()$summary
  }, digits = 6)

  output$exp_mean_dr <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(mean_dr ~ config_id + method, data = experiment_result()$summary))
  }, digits = 6)
  output$exp_mean_dr_norm <- shiny::renderTable({
    shiny::req(experiment_result())
    tab <- xtabs(mean_dr ~ config_id + method, data = experiment_result()$summary)
    matrix_to_df(normalize_xtab(tab))
  }, digits = 6)

  output$exp_median_dr <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(median_dr ~ config_id + method, data = experiment_result()$summary))
  }, digits = 6)
  output$exp_median_dr_norm <- shiny::renderTable({
    shiny::req(experiment_result())
    tab <- xtabs(median_dr ~ config_id + method, data = experiment_result()$summary)
    matrix_to_df(normalize_xtab(tab))
  }, digits = 6)

  output$exp_mean_dtheta <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(mean_dtheta ~ config_id + method, data = experiment_result()$summary))
  }, digits = 6)
  output$exp_mean_dtheta_norm <- shiny::renderTable({
    shiny::req(experiment_result())
    tab <- xtabs(mean_dtheta ~ config_id + method, data = experiment_result()$summary)
    matrix_to_df(normalize_xtab(tab))
  }, digits = 6)

  output$exp_median_dtheta <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(median_dtheta ~ config_id + method, data = experiment_result()$summary))
  }, digits = 6)
  output$exp_median_dtheta_norm <- shiny::renderTable({
    shiny::req(experiment_result())
    tab <- xtabs(median_dtheta ~ config_id + method, data = experiment_result()$summary)
    matrix_to_df(normalize_xtab(tab))
  }, digits = 6)

  output$exp_active_pixels <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(mean_active_pixels ~ config_id + method, data = experiment_result()$summary))
  }, digits = 2)
  output$exp_active_weight <- shiny::renderTable({
    shiny::req(experiment_result())
    matrix_to_df(xtabs(mean_active_weight ~ config_id + method, data = experiment_result()$summary))
  }, digits = 2)

  find_config <- shiny::reactive({
    if (input$find_config_source == "standard") {
      cfgs <- standard_single_configs(n_row = input$find_n_row, n_col = input$find_n_col)
      return(cfgs[[match(input$find_standard_config, vapply(cfgs, `[[`, character(1), "config_id"))]])
    }

    cfgs <- build_config_list(
      num_lines = input$find_num_lines,
      n_row = input$find_n_row,
      n_col = input$find_n_col,
      use_standard_single = FALSE,
      line_text = input$find_line_text
    )
    cfgs[[1]]
  })

  find_result <- shiny::eventReactive(input$run_find, {
    cfg <- find_config()
    method_map <- setNames(names(method_choices), unname(method_choices))
    method_name <- unname(method_map[input$find_method])
    method_name <- resolve_cssa_method_names(method_name, input$find_cssa_esprit)
    num_lines <- if (input$find_config_source == "standard") 1L else as.integer(input$find_num_lines)

    shiny::withProgress(message = "Поиск больших ошибок", value = 0, {
      out <- find_big_error_cases(
        cfg = cfg,
        method_name = method_name,
        sigma_noise = input$find_sigma,
        num_of_lines = num_lines,
        num_maxima = input$find_num_maxima,
        rho_step_ht = input$find_rho_step,
        theta_step_ht = input$find_theta_step,
        ht_type = input$find_ht_type,
        threshold_mode = input$find_threshold_mode,
        threshold_value = input$find_threshold_value,
        threshold_quantile_p = input$find_threshold_quantile_p,
        threshold_multiplier = input$find_threshold_multiplier,
        n_search = input$find_n_search,
        factor_threshold = input$find_factor,
        max_cases = input$find_max_cases,
        n_workers = input$find_n_workers,
        max_row_k = input$find_max_row_k
      )
      out
    })
  })

  output$find_ideal_table <- shiny::renderTable({
    shiny::req(find_result())
    ideal <- find_result()$ideal
    data.frame(
      ideal_dr = ideal$dr,
      ideal_dtheta = ideal$dtheta,
      baseline_dr = ideal$baseline_dr,
      baseline_dtheta = ideal$baseline_dtheta
    )
  }, digits = 6)

  output$find_cases_table <- shiny::renderTable({
    shiny::req(find_result())
    cases <- find_result()$cases
    if (length(cases) == 0L) {
      return(data.frame(message = "Случаи не найдены"))
    }
    dplyr::bind_rows(lapply(cases, function(case) {
      data.frame(
        case_id = case$case_id,
        rep = case$rep,
        dr = case$dr,
        dtheta = case$dtheta,
        dr_ratio = case$dr_ratio,
        dtheta_ratio = case$dtheta_ratio,
        threshold_value = case$threshold_value
      )
    }))
  }, digits = 6)

  output$find_case_selector <- shiny::renderUI({
    shiny::req(find_result())
    cases <- find_result()$cases
    if (length(cases) == 0L) {
      return(NULL)
    }
    shiny::selectInput(
      "find_case_id",
      "Случай для просмотра",
      choices = vapply(cases, `[[`, character(1), "case_id")
    )
  })

  selected_case <- shiny::reactive({
    shiny::req(find_result())
    cases <- find_result()$cases
    if (length(cases) == 0L) {
      return(NULL)
    }
    ids <- vapply(cases, `[[`, character(1), "case_id")
    cases[[match(input$find_case_id %||% ids[1], ids)]]
  })

  output$find_selected_case_meta <- shiny::renderTable({
    case <- selected_case()
    shiny::req(case)
    data.frame(
      case_id = case$case_id,
      rep = case$rep,
      dr = case$dr,
      dtheta = case$dtheta,
      dr_ratio = case$dr_ratio,
      dtheta_ratio = case$dtheta_ratio,
      threshold_value = case$threshold_value
    )
  }, digits = 6)

  output$find_processed_plot <- shiny::renderPlot({
    case <- selected_case()
    shiny::req(case)
    draw_case_matrix_triplet(case)
  })

  output$find_accumulator_plot <- shiny::renderPlot({
    case <- selected_case()
    shiny::req(case)
    draw_accumulator(case)
  })
}

shiny::shinyApp(ui = ui, server = server)
