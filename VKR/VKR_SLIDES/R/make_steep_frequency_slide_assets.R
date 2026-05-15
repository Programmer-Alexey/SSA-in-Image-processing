find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(cur, "SSA-in-Image-processing.Rproj"))) {
      return(cur)
    }
    next_dir <- dirname(cur)
    if (identical(next_dir, cur)) break
    cur <- next_dir
  }
  stop("РќРµ СѓРґР°Р»РѕСЃСЊ РЅР°Р№С‚Рё РєРѕСЂРµРЅСЊ РїСЂРѕРµРєС‚Р°.")
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
}
script_dir <- if (!is.null(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
project_root <- find_project_root(script_dir)
slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES")
tables_dir <- file.path(slides_dir, "tables")
images_dir <- file.path(slides_dir, "images")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

samples_path <- file.path(tables_dir, "steep_frequency_error_samples.csv")
if (!file.exists(samples_path)) {
  stop("РќРµС‚ С„Р°Р№Р»Р° ", samples_path, ". РЎРЅР°С‡Р°Р»Р° Р·Р°РїСѓСЃС‚РёС‚Рµ make_steep_frequency_error_table.R.")
}

samples <- read.csv(samples_path, stringsAsFactors = FALSE)
steep <- subset(samples, config_id == "steep_full")
steep$omega_normalized <- abs(steep$omega_wrapped) / (2 * pi)

first_rows <- aggregate(
  cbind(argmax_error, mse, near_mass_share, active_share_0001, active_share_01) ~ row_id + true_col + omega_normalized,
  data = steep,
  FUN = mean
)
first_rows <- first_rows[order(first_rows$row_id), ]
first_rows <- first_rows[first_rows$row_id %in% seq(1L, 19L, by = 2L), ]

slide_table <- data.frame(
  "$j$" = as.character(first_rows$true_col),
  "$\\omega_j$" = first_rows$omega_normalized,
  "$|\\hat j-j|$" = first_rows$argmax_error,
  check.names = FALSE
)

write_df_tex <- function(df, file, digits = 3L) {
  df_out <- df
  for (j in seq_along(df_out)) {
    if (is.numeric(df_out[[j]])) {
      df_out[[j]] <- format(round(df_out[[j]], digits), nsmall = digits, trim = TRUE)
    } else {
      df_out[[j]] <- as.character(df_out[[j]])
      df_out[[j]] <- gsub("_", "\\\\_", df_out[[j]])
    }
  }
  align <- paste(rep("l", ncol(df_out)), collapse = "|")
  lines <- c(
    paste0("\\begin{tabular}{|", align, "|}"),
    "\\hline",
    paste0(paste(names(df_out), collapse = " & "), " \\\\"),
    "\\hline"
  )
  for (i in seq_len(nrow(df_out))) {
    lines <- c(lines, paste0(paste(df_out[i, ], collapse = " & "), " \\\\"), "\\hline")
  }
  lines <- c(lines, "\\end{tabular}")
  writeLines(lines, file, useBytes = TRUE)
}

write_df_tex(slide_table, file.path(tables_dir, "steep_full_first_rows_frequency.tex"), digits = 3L)

by_frequency <- aggregate(
  cbind(argmax_error, mse) ~ true_col + omega_normalized,
  data = steep,
  FUN = mean
)
by_frequency <- by_frequency[order(by_frequency$omega_normalized), ]
write.csv(by_frequency, file.path(tables_dir, "steep_full_frequency_by_row.csv"), row.names = FALSE)

pdf(file.path(images_dir, "steep_full_frequency_localization.pdf"), width = 5.3, height = 3.5)
plot(
  by_frequency$omega_normalized,
  by_frequency$argmax_error,
  type = "b",
  pch = 19,
  col = "#1f4e79",
  xlab = expression(omega[j] == (j - 1) / N),
  ylab = expression(abs(hat(j) - j)),
  main = "Mean localization error"
)
grid(col = "gray85")
abline(lm(argmax_error ~ omega_normalized, data = by_frequency), col = "red", lwd = 2)
dev.off()

png(file.path(images_dir, "steep_full_frequency_localization.png"), width = 1100, height = 760, res = 180)
plot(
  by_frequency$omega_normalized,
  by_frequency$argmax_error,
  type = "b",
  pch = 19,
  col = "#1f4e79",
  xlab = expression(omega[j] == (j - 1) / N),
  ylab = expression(abs(hat(j) - j)),
  main = "Mean localization error"
)
grid(col = "gray85")
abline(lm(argmax_error ~ omega_normalized, data = by_frequency), col = "red", lwd = 2)
dev.off()

cat("Р“РѕС‚РѕРІРѕ: steep_full_first_rows_frequency.tex Рё steep_full_frequency_localization.pdf СЃ РЅРѕСЂРјРёСЂРѕРІР°РЅРЅРѕР№ С‡Р°СЃС‚РѕС‚РѕР№\n")

