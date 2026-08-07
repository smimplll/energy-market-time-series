# Energy market time-series analysis
#
# The script runs two connected studies in sequence:
#   1. Univariate ARIMA analysis of Brent log returns.
#   2. VAR analysis of Brent, Henry Hub gas and the broad U.S. dollar index.
#
# Run from the repository root with:
#   Rscript R/time_series_analysis.R

options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA_RAW_DIR <- file.path(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures")
TABLE_DIR <- file.path(PROJECT_ROOT, "results", "tables")

for (path in c(DATA_RAW_DIR, DATA_PROCESSED_DIR, FIGURE_DIR, TABLE_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("forecast", "ggplot2", "tseries", "urca", "vars")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required packages first with source(\"install_packages.R\"): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

START_DATE <- as.Date("2010-01-04")
END_DATE <- as.Date("2025-11-24")
REFRESH_DATA <- tolower(Sys.getenv("REFRESH_DATA", "true")) %in%
  c("1", "true", "yes", "y")

fred_url <- function(series_id, start_date = START_DATE, end_date = END_DATE) {
  sprintf(
    paste0(
      "https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s",
      "&cosd=%s&coed=%s"
    ),
    series_id,
    format(start_date, "%Y-%m-%d"),
    format(end_date, "%Y-%m-%d")
  )
}

read_fred_series <- function(
  series_id,
  value_name,
  start_date = START_DATE,
  end_date = END_DATE,
  refresh = REFRESH_DATA
) {
  cache_file <- file.path(DATA_RAW_DIR, paste0(series_id, ".csv"))
  source_url <- fred_url(series_id, start_date, end_date)

  if (refresh) {
    temporary_file <- tempfile(fileext = ".csv")
    download_ok <- tryCatch(
      {
        status <- utils::download.file(
          source_url,
          temporary_file,
          mode = "wb",
          quiet = TRUE
        )
        identical(status, 0L) && file.info(temporary_file)$size > 0
      },
      error = function(error) {
        warning(
          "FRED download failed for ", series_id,
          "; the cached file will be used if available. Details: ",
          conditionMessage(error),
          call. = FALSE
        )
        FALSE
      }
    )

    if (isTRUE(download_ok)) {
      file.copy(temporary_file, cache_file, overwrite = TRUE)
    }
  }

  if (!file.exists(cache_file)) {
    stop(
      "No downloaded or cached data found for ", series_id,
      ". Check the internet connection or add ", cache_file,
      call. = FALSE
    )
  }

  data <- utils::read.csv(
    cache_file,
    na.strings = c("", ".", "NA", "N/A"),
    check.names = FALSE
  )

  if (!all(c("observation_date", series_id) %in% names(data))) {
    stop("Unexpected FRED file structure for ", series_id, call. = FALSE)
  }

  result <- data.frame(
    date = as.Date(data$observation_date),
    value = suppressWarnings(as.numeric(data[[series_id]]))
  )
  names(result)[2] <- value_name

  result <- result[
    !is.na(result$date) & result$date >= start_date & result$date <= end_date,
    ,
    drop = FALSE
  ]
  result <- result[order(result$date), , drop = FALSE]
  attr(result, "source_url") <- source_url
  result
}

with_png <- function(filename, plot_code, width = 1800, height = 1100, res = 180) {
  path <- file.path(FIGURE_DIR, filename)
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(plot_code)
  invisible(path)
}

save_ggplot <- function(plot, filename, width = 9, height = 5.5) {
  path <- file.path(FIGURE_DIR, filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 180)
  invisible(path)
}

capture_to_file <- function(object, filename) {
  capture.output(object, file = file.path(TABLE_DIR, filename))
}

message("Downloading FRED series or reading cached snapshots...")
brent_raw <- read_fred_series("DCOILBRENTEU", "brent")
gas_raw <- read_fred_series("DHHNGSP", "gas")
dollar_raw <- read_fred_series("DTWEXBGS", "dollar_index")

source_manifest <- data.frame(
  series_id = c("DCOILBRENTEU", "DHHNGSP", "DTWEXBGS"),
  variable = c("brent", "gas", "dollar_index"),
  source_url = c(
    attr(brent_raw, "source_url"),
    attr(gas_raw, "source_url"),
    attr(dollar_raw, "source_url")
  ),
  start_date = START_DATE,
  end_date = END_DATE
)
utils::write.csv(
  source_manifest,
  file.path(TABLE_DIR, "00_data_sources.csv"),
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------
# Part 1. Univariate ARIMA analysis of Brent
# -----------------------------------------------------------------------------

message("Part 1: ARIMA analysis of Brent...")

brent <- brent_raw[!is.na(brent_raw$brent) & brent_raw$brent > 0, , drop = FALSE]
brent$log_price <- log(brent$brent)
brent$log_return <- c(NA_real_, diff(brent$log_price))
brent_returns <- stats::na.omit(brent$log_return)
brent_return_dates <- brent$date[!is.na(brent$log_return)]

brent_ts <- stats::ts(brent$brent, frequency = 252)
brent_log_ts <- stats::ts(brent$log_price, frequency = 252)
brent_return_ts <- stats::ts(brent_returns, frequency = 252)

brent_summary <- data.frame(
  raw_observations = nrow(brent_raw),
  missing_prices = sum(is.na(brent_raw$brent)),
  usable_prices = nrow(brent),
  first_date = min(brent$date),
  last_date = max(brent$date),
  mean_price = mean(brent$brent),
  standard_deviation_price = stats::sd(brent$brent),
  mean_log_return = mean(brent_returns),
  standard_deviation_log_return = stats::sd(brent_returns)
)
utils::write.csv(
  brent_summary,
  file.path(TABLE_DIR, "01_brent_summary.csv"),
  row.names = FALSE
)

price_plot <- ggplot2::ggplot(brent, ggplot2::aes(x = date, y = brent)) +
  ggplot2::geom_line(color = "#9C2C2C", linewidth = 0.45) +
  ggplot2::labs(
    title = "Brent spot price",
    subtitle = "Europe Brent Spot Price FOB",
    x = NULL,
    y = "USD per barrel"
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_ggplot(price_plot, "01_brent_price.png")

return_plot <- ggplot2::ggplot(
  data.frame(date = brent_return_dates, log_return = brent_returns),
  ggplot2::aes(x = date, y = log_return)
) +
  ggplot2::geom_line(color = "#285F8F", linewidth = 0.35) +
  ggplot2::labs(
    title = "Brent logarithmic returns",
    x = NULL,
    y = "Log return"
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_ggplot(return_plot, "02_brent_log_returns.png")

with_png(
  "03_brent_price_acf_pacf.png",
  forecast::tsdisplay(brent_ts, main = "Brent price: series, ACF and PACF")
)
with_png(
  "04_brent_log_price_acf_pacf.png",
  forecast::tsdisplay(brent_log_ts, main = "Log Brent price: series, ACF and PACF")
)
with_png(
  "05_brent_returns_acf_pacf.png",
  forecast::tsdisplay(brent_return_ts, main = "Brent log returns: series, ACF and PACF")
)

adf_price <- urca::ur.df(
  brent_log_ts,
  type = "trend",
  lags = 10,
  selectlags = "AIC"
)
kpss_price <- tseries::kpss.test(brent_log_ts)
adf_returns <- urca::ur.df(
  brent_return_ts,
  type = "drift",
  lags = 10,
  selectlags = "AIC"
)
kpss_returns <- tseries::kpss.test(brent_return_ts)

capture_to_file(summary(adf_price), "02_adf_log_price.txt")
capture_to_file(kpss_price, "03_kpss_log_price.txt")
capture_to_file(summary(adf_returns), "04_adf_log_returns.txt")
capture_to_file(kpss_returns, "05_kpss_log_returns.txt")

candidate_orders <- rbind(
  c(0, 0, 0), c(1, 0, 0), c(2, 0, 0), c(0, 0, 1),
  c(0, 0, 2), c(1, 0, 1), c(1, 0, 2), c(2, 0, 1)
)

candidate_models <- lapply(seq_len(nrow(candidate_orders)), function(i) {
  forecast::Arima(
    brent_return_ts,
    order = candidate_orders[i, ],
    include.mean = TRUE,
    method = "ML"
  )
})

arima_comparison <- data.frame(
  model = apply(
    candidate_orders,
    1,
    function(order) sprintf("ARIMA(%d,%d,%d)", order[1], order[2], order[3])
  ),
  AIC = vapply(candidate_models, stats::AIC, numeric(1)),
  BIC = vapply(candidate_models, stats::BIC, numeric(1))
)

zero_mean_model <- forecast::Arima(
  brent_return_ts,
  order = c(0, 0, 0),
  include.mean = FALSE,
  method = "ML"
)
arima_comparison <- rbind(
  arima_comparison,
  data.frame(
    model = "ARIMA(0,0,0), zero mean",
    AIC = stats::AIC(zero_mean_model),
    BIC = stats::BIC(zero_mean_model)
  )
)
arima_comparison <- arima_comparison[order(arima_comparison$BIC), ]
utils::write.csv(
  arima_comparison,
  file.path(TABLE_DIR, "06_arima_model_comparison.csv"),
  row.names = FALSE
)

automatic_aicc <- forecast::auto.arima(brent_return_ts)
automatic_bic <- forecast::auto.arima(brent_return_ts, ic = "bic")
final_arima <- zero_mean_model

capture_to_file(summary(automatic_aicc), "07_auto_arima_aicc.txt")
capture_to_file(summary(automatic_bic), "08_auto_arima_bic.txt")
capture_to_file(summary(final_arima), "09_final_arima.txt")

with_png(
  "06_arima_residual_diagnostics.png",
  forecast::checkresiduals(final_arima),
  width = 1800,
  height = 1400
)

test_horizon <- 250L
total_observations <- length(brent_return_ts)
if (total_observations <= test_horizon) {
  stop("The Brent return series is too short for the requested test horizon.")
}

training_indices <- seq_len(total_observations - test_horizon)
test_indices <- (total_observations - test_horizon + 1L):total_observations
brent_train <- brent_return_ts[training_indices]
brent_test <- brent_return_ts[test_indices]

training_arima <- forecast::Arima(
  brent_train,
  order = c(0, 0, 0),
  include.mean = FALSE,
  method = "ML"
)
brent_forecast <- forecast::forecast(training_arima, h = length(brent_test))
brent_accuracy <- forecast::accuracy(brent_forecast, brent_test)

with_png(
  "07_arima_out_of_sample_forecast.png",
  {
    plot(
      seq_along(brent_test),
      as.numeric(brent_forecast$mean),
      type = "l",
      col = "#285F8F",
      lwd = 2,
      xlab = "Test observation",
      ylab = "Log return",
      main = "Brent log-return forecast: train/test",
      ylim = range(
        c(
          brent_forecast$lower[, 2],
          brent_forecast$upper[, 2],
          brent_test
        ),
        na.rm = TRUE
      )
    )
    lines(brent_forecast$lower[, 2], col = "#285F8F", lty = 2)
    lines(brent_forecast$upper[, 2], col = "#285F8F", lty = 2)
    lines(as.numeric(brent_test), col = "#B83A3A", lwd = 1)
    grid()
    legend(
      "topright",
      legend = c("Forecast", "Observed"),
      col = c("#285F8F", "#B83A3A"),
      lty = 1,
      lwd = c(2, 1),
      bty = "n"
    )
  }
)

set.seed(123)
innovation_sd <- sqrt(final_arima$sigma2)
synthetic_returns <- stats::ts(
  stats::rnorm(length(brent_return_ts), mean = 0, sd = innovation_sd),
  frequency = 252
)

adf_synthetic <- urca::ur.df(
  synthetic_returns,
  type = "drift",
  selectlags = "AIC"
)
kpss_synthetic <- tseries::kpss.test(synthetic_returns)
automatic_synthetic <- forecast::auto.arima(synthetic_returns)

capture_to_file(summary(adf_synthetic), "10_adf_synthetic.txt")
capture_to_file(kpss_synthetic, "11_kpss_synthetic.txt")
capture_to_file(summary(automatic_synthetic), "12_auto_arima_synthetic.txt")

with_png(
  "08_synthetic_returns_acf_pacf.png",
  forecast::tsdisplay(
    synthetic_returns,
    main = "Synthetic Brent log returns: series, ACF and PACF"
  )
)

synthetic_train <- synthetic_returns[training_indices]
synthetic_test <- synthetic_returns[test_indices]
synthetic_training_model <- forecast::Arima(
  synthetic_train,
  order = c(0, 0, 0),
  include.mean = FALSE,
  method = "ML"
)
synthetic_forecast <- forecast::forecast(
  synthetic_training_model,
  h = length(synthetic_test)
)
synthetic_accuracy <- forecast::accuracy(synthetic_forecast, synthetic_test)

with_png(
  "09_synthetic_out_of_sample_forecast.png",
  {
    plot(
      seq_along(synthetic_test),
      as.numeric(synthetic_forecast$mean),
      type = "l",
      col = "#285F8F",
      lwd = 2,
      xlab = "Test observation",
      ylab = "Synthetic log return",
      main = "Synthetic-series forecast: train/test",
      ylim = range(
        c(
          synthetic_forecast$lower[, 2],
          synthetic_forecast$upper[, 2],
          synthetic_test
        ),
        na.rm = TRUE
      )
    )
    lines(synthetic_forecast$lower[, 2], col = "#285F8F", lty = 2)
    lines(synthetic_forecast$upper[, 2], col = "#285F8F", lty = 2)
    lines(as.numeric(synthetic_test), col = "#B83A3A", lwd = 1)
    grid()
    legend(
      "topright",
      legend = c("Forecast", "Observed"),
      col = c("#285F8F", "#B83A3A"),
      lty = 1,
      lwd = c(2, 1),
      bty = "n"
    )
  }
)

accuracy_rows <- function(matrix, series_name) {
  data.frame(
    series = series_name,
    sample = rownames(matrix),
    RMSE = matrix[, "RMSE"],
    MAE = matrix[, "MAE"],
    row.names = NULL
  )
}

forecast_accuracy <- rbind(
  accuracy_rows(brent_accuracy, "Observed Brent returns"),
  accuracy_rows(synthetic_accuracy, "Synthetic Brent returns")
)
utils::write.csv(
  forecast_accuracy,
  file.path(TABLE_DIR, "13_arima_forecast_accuracy.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Part 2. Multivariate VAR analysis
# -----------------------------------------------------------------------------

message("Part 2: VAR analysis of Brent, Henry Hub gas and the dollar index...")

market_levels <- Reduce(
  function(left, right) merge(left, right, by = "date", all = FALSE),
  list(brent_raw, gas_raw, dollar_raw)
)
market_levels <- market_levels[stats::complete.cases(market_levels), ]
market_levels <- market_levels[
  market_levels$brent > 0 &
    market_levels$gas > 0 &
    market_levels$dollar_index > 0,
  ,
  drop = FALSE
]
market_levels <- market_levels[order(market_levels$date), ]

market_returns <- data.frame(
  date = market_levels$date[-1],
  r_brent = diff(log(market_levels$brent)),
  r_gas = diff(log(market_levels$gas)),
  r_dollar = diff(log(market_levels$dollar_index))
)
market_returns <- market_returns[stats::complete.cases(market_returns), ]

utils::write.csv(
  market_levels,
  file.path(DATA_PROCESSED_DIR, "market_levels.csv"),
  row.names = FALSE
)
utils::write.csv(
  market_returns,
  file.path(DATA_PROCESSED_DIR, "market_log_returns.csv"),
  row.names = FALSE
)

level_summary <- data.frame(
  observations = nrow(market_levels),
  first_date = min(market_levels$date),
  last_date = max(market_levels$date),
  mean_brent = mean(market_levels$brent),
  mean_gas = mean(market_levels$gas),
  mean_dollar_index = mean(market_levels$dollar_index)
)
utils::write.csv(
  level_summary,
  file.path(TABLE_DIR, "14_var_data_summary.csv"),
  row.names = FALSE
)

series_specs <- list(
  brent = list(color = "#397BB7", label = "Brent", unit = "USD per barrel"),
  gas = list(color = "#2C7A3F", label = "Henry Hub gas", unit = "USD per MMBtu"),
  dollar_index = list(color = "#8E44D6", label = "Broad U.S. dollar index", unit = "Index")
)

for (variable in names(series_specs)) {
  spec <- series_specs[[variable]]
  plot <- ggplot2::ggplot(
    market_levels,
    ggplot2::aes(x = date, y = .data[[variable]])
  ) +
    ggplot2::geom_line(color = spec$color, linewidth = 0.4) +
    ggplot2::labs(title = spec$label, x = NULL, y = spec$unit) +
    ggplot2::theme_minimal(base_size = 11)
  save_ggplot(plot, paste0("10_level_", variable, ".png"))
}

return_specs <- list(
  r_brent = list(color = "#397BB7", label = "Brent log returns"),
  r_gas = list(color = "#2C7A3F", label = "Henry Hub log returns"),
  r_dollar = list(color = "#8E44D6", label = "Dollar-index log returns")
)

for (variable in names(return_specs)) {
  spec <- return_specs[[variable]]
  plot <- ggplot2::ggplot(
    market_returns,
    ggplot2::aes(x = date, y = .data[[variable]])
  ) +
    ggplot2::geom_line(color = spec$color, linewidth = 0.3) +
    ggplot2::labs(title = spec$label, x = NULL, y = "Log return") +
    ggplot2::theme_minimal(base_size = 11)
  save_ggplot(plot, paste0("11_return_", variable, ".png"))
}

stationarity_tests <- lapply(
  market_returns[c("r_brent", "r_gas", "r_dollar")],
  function(series) list(
    adf = tseries::adf.test(series),
    kpss = tseries::kpss.test(series)
  )
)
capture_to_file(stationarity_tests, "15_var_stationarity_tests.txt")

return_ts <- stats::ts(
  market_returns[c("r_brent", "r_gas", "r_dollar")],
  frequency = 252
)

lag_selection <- vars::VARselect(return_ts, lag.max = 10, type = "const")
capture_to_file(lag_selection, "16_var_lag_selection.txt")

selected_lags <- 2L
var_model <- vars::VAR(return_ts, p = selected_lags, type = "const")
capture_to_file(summary(var_model), "17_var_model_summary.txt")

var_roots <- vars::roots(var_model, modulus = TRUE)
utils::write.csv(
  data.frame(root = seq_along(var_roots), modulus = var_roots),
  file.path(TABLE_DIR, "18_var_root_moduli.csv"),
  row.names = FALSE
)

with_png(
  "12_var_recursive_cusum.png",
  plot(vars::stability(var_model, type = "Rec-CUSUM")),
  width = 1800,
  height = 1500
)

serial_diagnostic <- vars::serial.test(
  var_model,
  lags.pt = 16,
  type = "PT.asymptotic"
)
arch_diagnostic <- vars::arch.test(var_model, lags.multi = 5)
normality_diagnostic <- vars::normality.test(var_model)

capture_to_file(serial_diagnostic, "19_var_serial_test.txt")
capture_to_file(arch_diagnostic, "20_var_arch_test.txt")
capture_to_file(normality_diagnostic, "21_var_normality_test.txt")

granger_results <- list(
  dollar = vars::causality(var_model, cause = "r_dollar"),
  gas = vars::causality(var_model, cause = "r_gas"),
  brent = vars::causality(var_model, cause = "r_brent")
)
capture_to_file(granger_results, "22_var_granger_causality.txt")

set.seed(123)
irf_dollar_to_brent <- vars::irf(
  var_model,
  impulse = "r_dollar",
  response = "r_brent",
  n.ahead = 20,
  ortho = TRUE,
  boot = TRUE,
  runs = 500,
  ci = 0.95
)
irf_gas_to_brent <- vars::irf(
  var_model,
  impulse = "r_gas",
  response = "r_brent",
  n.ahead = 20,
  ortho = TRUE,
  boot = TRUE,
  runs = 500,
  ci = 0.95
)
irf_brent_to_gas <- vars::irf(
  var_model,
  impulse = "r_brent",
  response = "r_gas",
  n.ahead = 20,
  ortho = TRUE,
  boot = TRUE,
  runs = 500,
  ci = 0.95
)

with_png("13_irf_dollar_to_brent.png", plot(irf_dollar_to_brent))
with_png("14_irf_gas_to_brent.png", plot(irf_gas_to_brent))
with_png("15_irf_brent_to_gas.png", plot(irf_brent_to_gas))

variance_decomposition <- vars::fevd(var_model, n.ahead = 12)
with_png(
  "16_var_forecast_error_variance_decomposition.png",
  plot(variance_decomposition),
  width = 1800,
  height = 1500
)

var_forecast <- stats::predict(var_model, n.ahead = 12, ci = 0.95)
with_png(
  "17_var_forecast.png",
  plot(var_forecast, xlab = "Forecast horizon"),
  width = 1800,
  height = 1500
)

brent_var_forecast <- as.data.frame(var_forecast$fcst$r_brent)
brent_var_forecast$horizon <- seq_len(nrow(brent_var_forecast))
brent_var_forecast <- brent_var_forecast[
  c("horizon", setdiff(names(brent_var_forecast), "horizon"))
]
utils::write.csv(
  brent_var_forecast,
  file.path(TABLE_DIR, "23_var_brent_forecast.csv"),
  row.names = FALSE
)

capture.output(
  utils::sessionInfo(),
  file = file.path(TABLE_DIR, "24_session_info.txt")
)

message("Analysis completed. Results are in results/figures and results/tables.")
