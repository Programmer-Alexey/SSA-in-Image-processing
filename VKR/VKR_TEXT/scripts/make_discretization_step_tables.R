script_dir_for_source <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source(file.path(script_dir_for_source(), "00_common.R"))

ensure_report_dirs()
load_project_code(load_hough = TRUE)

n_row <- 100L
n_col <- 100L
line <- data.frame(a = 2, b = -1, intensity = 1)
clean <- make_line_image(n_row, n_col, line, line_method = "default")
true_line <- matrix(convert_ab_to_rho_theta(line$a, line$b), nrow = 1L)
colnames(true_line) <- c("rho", "theta")

grid_df <- data.frame(
  rho_step = c(5, 1, 0.5, 0.25, 0.05, 0.02),
  theta_step = c(1, 0.01, 0.005, 0.0025, 0.005, 0.002),
  stringsAsFactors = FALSE
)

samples <- do.call(rbind, lapply(seq_len(nrow(grid_df)), function(i) {
  pred <- project_to_hough_grid(
    true_line,
    rho_step = grid_df$rho_step[i],
    theta_step = grid_df$theta_step[i]
  )
  err <- parameter_line_errors_by_line(true_line, pred)
  data.frame(
    rho_step = grid_df$rho_step[i],
    theta_step = grid_df$theta_step[i],
    pred_rho = pred[1L, "rho"],
    pred_theta = pred[1L, "theta"],
    rho_mse = err$rho_mse,
    theta_mse = err$theta_mse,
    stringsAsFactors = FALSE
  )
}))

write_csv_utf8(samples, file.path(chapter3_data_dir, "discretization_step_summary.csv"))

tex_df <- data.frame(
  "$h_\\rho$" = format_plain_num(samples$rho_step, 4L),
  "$h_\\theta$" = format_plain_num(samples$theta_step, 4L),
  "MSE по $\\rho$" = ifelse(samples$rho_mse < 1e-3, sprintf("%.2e", samples$rho_mse), format_plain_num(samples$rho_mse, 6L)),
  "MSE по $\\theta$" = ifelse(samples$theta_mse < 1e-3, sprintf("%.2e", samples$theta_mse), format_plain_num(samples$theta_mse, 6L)),
  check.names = FALSE
)

write_lines_utf8(
  latex_table_lines(tex_df, resize = NULL),
  file.path(chapter3_table_dir, "chapter3_discretization_decay.tex")
)

cat("Done. Discretization table saved to: ",
    file.path(chapter3_table_dir, "chapter3_discretization_decay.tex"),
    "\n",
    sep = "")
