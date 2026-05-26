script_dir_for_source <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- script_dir_for_source()
source(file.path(script_dir, "make_report_figures.R"))

