packages <- c("forecast", "ggplot2", "tseries", "urca", "vars")

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, dependencies = TRUE)
}

message("All required R packages are installed.")
