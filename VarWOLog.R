# =========================
# Libraries & Helpers
# =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(forecast)
  library(tseries)
  library(ggplot2)
  library(vars)
  library(lmtest)
})

safe_accuracy <- function(fc, ytest) {
  A <- accuracy(fc, ytest)
  A[, intersect(c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1","Theil's U"),
                colnames(A)), drop = FALSE]
}

lb_p <- function(fit, y, lag = 8) {
  k <- length(coef(fit)); if (is.null(k)) k <- 0
  Box.test(residuals(fit), type = "Ljung-Box",
           lag = min(lag, length(y)-1), fitdf = k)$p.value
}

get_AICc <- function(fit) {  # robust AICc
  if (!is.null(fit$aicc)) return(fit$aicc)
  aic <- AIC(fit); n <- length(residuals(fit)); k <- length(coef(fit))
  if (is.finite(n) && is.finite(k) && (n > k + 1))
    aic + (2 * k * (k + 1)) / (n - k - 1) else aic
}

plot_backtest_line <- function(ts_full, yhat_test, start_yq, title, ylab = NULL) {
  n  <- length(ts_full)
  f  <- frequency(ts_full)
  st <- start(ts_full)
  all_t <- st[1] + (0:(n - 1) + st[2] - 1) / f
  
  q_to_month <- function(q) c(1, 4, 7, 10)[q]
  yy <- floor(all_t)
  qq <- (round((all_t - yy) * f) %% f) + 1
  dates_full <- as.Date(sprintf("%04d-%02d-01", yy, q_to_month(qq)))
  
  test_actual <- window(ts_full, start = start_yq)
  idx_test    <- (n - length(test_actual) + 1):n
  
  if (length(test_actual) != length(yhat_test)) {
    stop("Lengths of actual test series and forecasts do not match.")
  }
  
  df <- tibble(
    date   = dates_full[idx_test],
    actual = as.numeric(test_actual),
    mean   = as.numeric(yhat_test)
  )
  
  ggplot(df, aes(date, actual)) +
    geom_line(linewidth = 0.9) +
    geom_line(aes(y = mean), linewidth = 0.9, linetype = "dashed") +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 12)
}

adf_k <- function(x) trunc((length(na.omit(x)) - 1)^(1/3))


# =========================
# Load & Clean Data
# =========================
df <- read.csv("DAta.csv", stringsAsFactors = FALSE)

norm_str <- function(x) tolower(trimws(as.character(x)))
df <- df %>%
  mutate(
    series = norm_str(series),
    sector = norm_str(sector),
    value  = as.numeric(gsub(",", "", as.character(value)))
  )

d0 <- suppressWarnings(dmy(df$date))
d1 <- suppressWarnings(ymd(df$date))
d2 <- suppressWarnings(mdy(df$date))
df$date <- coalesce(d0, d1, d2)
stopifnot(!all(is.na(df$date)))

abs_df <- df %>%
  filter(series == "abs", !is.na(value)) %>%
  arrange(date)

# =========================
# Build ts_p0, ts_p1, ts_share
# =========================
gdp  <- abs_df %>%
  dplyr::filter(sector %in% c("p0","gdp","gross domestic product","jumlah","kdnk")) %>%
  dplyr::arrange(date) %>%
  dplyr::select(date, p0 = value)

agri <- abs_df %>%
  dplyr::filter(sector %in% c("p1","agriculture","agriculture, forestry & fishing","pertanian")) %>%
  dplyr::arrange(date) %>%
  dplyr::select(date, p1 = value)

dat <- inner_join(gdp, agri, by = "date") %>%
  arrange(date)

ts_p0 <- ts(dat$p0,
            start = c(year(min(dat$date)), quarter(min(dat$date))),
            frequency = 4)
ts_p1 <- ts(dat$p1, start = start(ts_p0), frequency = 4)
ts_share <- ts(100 * dat$p1 / dat$p0,
               start = start(ts_p0), frequency = 4)

# Optional plots
plot_dat <- dat %>%
  mutate(share = 100 * p1/p0) %>%
  rename(`GDP (p0)` = p0,
         `Agriculture (p1)` = p1,
         `Share of GDP (%)` = share)

p_facets <- plot_dat |>
  pivot_longer(cols = c(`GDP (p0)`, `Agriculture (p1)`, `Share of GDP (%)`),
               names_to = "series", values_to = "value") |>
  ggplot(aes(date, value)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  labs(title = "GDP, Agriculture, and Share over time", x = NULL, y = NULL) +
  theme_minimal(base_size = 12)

# print(p_facets)


# =========================
# Pulses / Intervention Dummies (xreg)
# =========================
tq <- time(ts_p0)

X <- cbind(
  covid_q2     = as.numeric(tq == 2020 + 1/4),  # 2020Q2
  covid_q3     = as.numeric(tq == 2020 + 2/4),  # 2020Q3
  shock_2022q2 = as.numeric(tq == 2022 + 1/4)   # 2022Q2
)

X_ts <- ts(X, start = start(ts_p0), frequency = 4)


# =========================
# Stationarity Checks (ADF)
# =========================
log_p0 <- log(ts_p0)
log_p1 <- log(ts_p1)
log_share <- log(ts_p1) - log(ts_p0)  # log(p1/p0)

# Raw levels (no diff)
adf_p0_lvl    <- adf.test(ts_p0,    k = adf_k(ts_p0))
adf_p1_lvl    <- adf.test(ts_p1,    k = adf_k(ts_p1))
adf_share_lvl <- adf.test(ts_share, k = adf_k(ts_share))

# Raw: seasonal differences (lag 4)
adf_p0_sdiff    <- adf.test(diff(ts_p0,    lag = 4), k = adf_k(ts_p0))
adf_p1_sdiff    <- adf.test(diff(ts_p1,    lag = 4), k = adf_k(ts_p1))
adf_share_sdiff <- adf.test(diff(ts_share, lag = 4), k = adf_k(ts_share))

# Logs: levels
adf_logp0_lvl    <- adf.test(log_p0,    k = adf_k(log_p0))
adf_logp1_lvl    <- adf.test(log_p1,    k = adf_k(log_p1))
adf_logshare_lvl <- adf.test(log_share, k = adf_k(log_share))

# Logs: seasonal differences (lag 4)
adf_logp0_sdiff    <- adf.test(diff(log_p0,    lag = 4), k = adf_k(log_p0))
adf_logp1_sdiff    <- adf.test(diff(log_p1,    lag = 4), k = adf_k(log_p1))
adf_logshare_sdiff <- adf.test(diff(log_share, lag = 4), k = adf_k(log_share))

# (You can print these if you like)
# adf_logp0_sdiff; adf_logp1_sdiff; adf_logshare_sdiff


# =========================
# Correlation p0 vs p1
# =========================
cor_p0_p1_full <- cor(as.numeric(ts_p0), as.numeric(ts_p1), use = "complete.obs")

train_p0 <- window(ts_p0, end = c(2023, 4))
train_p1 <- window(ts_p1, end = c(2023, 4))
cor_p0_p1_train <- cor(as.numeric(train_p0), as.numeric(train_p1), use = "complete.obs")

d4_p0 <- diff(ts_p0, lag = 4)
d4_p1 <- diff(ts_p1, lag = 4)
cor_d4_p0_p1 <- cor(as.numeric(d4_p0), as.numeric(d4_p1), use = "complete.obs")

# =========================
# VAR on Seasonal Differences of Logs + xreg
# =========================

# 1) Build log models
d4_log_p0 <- diff(log_p0, lag = 4)  # Δ4 log(p0)
d4_log_p1 <- diff(log_p1, lag = 4)  # Δ4 log(p1) not used currently

y_var <- ts(
  cbind(log_p0, log_p1),
  start     = start(log_p0),
  frequency = 4
)
colnames(y_var) <- c("log_p0", "log_p1")

# Align xreg pulses with y_var
X_var <- window(X_ts, start = start(y_var), end = end(y_var))
stopifnot(nrow(X_var) == nrow(y_var))

# 2) Train–test split (same as ARIMA): train ≤ 2023Q4, test = 2024Q1–2025Q4
cut_year <- c(2023, 4)

y_train <- window(y_var, end   = cut_year)
y_test  <- window(y_var, start = c(2024, 1))

X_train <- window(X_var, end   = cut_year)
X_test  <- window(X_var, start = c(2024, 1))

stopifnot(
  nrow(y_train) == nrow(X_train),
  nrow(y_test)  == nrow(X_test)
)

# 3) Lag selection on TRAIN only
lag_sel <- VARselect(y_train, lag.max = 4, type = "trend", exogen = X_train)
lag_sel



p_opt <- lag_sel$selection["AIC(n)"]  # you can switch to "SC(n)" if you prefer BIC

jo_tr <- ca.jo(y_train, type = "trace", K = p_opt, ecdet = "trend")
summary(jo_tr)

# 4) Fit VAR on TRAIN (Δ4 log + xreg)
var_fit <- VAR(
  y_train,
  p      = p_opt,
  type   = "const",
  exogen = X_train
)
summary(var_fit)

# Diagnostics
serial.test(var_fit, lags.pt = 16, type = "PT.asymptotic")
normality.test(var_fit)
arch.test(var_fit, lags.multi = 4)
roots(var_fit)

# 5) Forecast over TEST horizon (2024–2025)
h <- nrow(y_test)

var_fc <- predict(
  var_fit,
  n.ahead = h,
  dumvar  = X_test
)

fc_d4_log_p0 <- ts(var_fc$fcst$log_p0[, "fcst"],
                   start = start(y_test), frequency = 4)
fc_d4_log_p1 <- ts(var_fc$fcst$log_p1[, "fcst"],
                   start = start(y_test), frequency = 4)

# True Δ4 logs in test period
test_d4_log_p0 <- window(log_p0, start = c(2024, 1))
test_d4_log_p1 <- window(log_p1, start = c(2024, 1))

cat("\n== TEST ACCURACY – VAR(Δ4 log + xreg) on Δ4 log scale ==\n")
cat("\nΔ4 log p0:\n")
print(safe_accuracy(fc_d4_log_p0, test_d4_log_p0))
cat("\nΔ4 log p1:\n")
print(safe_accuracy(fc_d4_log_p1, test_d4_log_p1))


# =========================
# Reconstruct Log Levels and Levels from Δ4 Log Forecasts
# =========================

# Anchors: last 4 observed logs up to 2023Q4 (end of train)
last_4_log_p0 <- tail(log_p0[time(log_p0) <= 2023 + 3/4], 4)
last_4_log_p1 <- tail(log_p1[time(log_p1) <= 2023 + 3/4], 4)

h <- length(fc_d4_log_p0)

fc_log_p0_full <- numeric(h + 4)
fc_log_p1_full <- numeric(h + 4)

fc_log_p0_full[1:4] <- as.numeric(last_4_log_p0)
fc_log_p1_full[1:4] <- as.numeric(last_4_log_p1)

for (k in 1:h) {
  fc_log_p0_full[k + 4] <- as.numeric(fc_d4_log_p0[k]) + fc_log_p0_full[k]
  fc_log_p1_full[k + 4] <- as.numeric(fc_d4_log_p1[k]) + fc_log_p1_full[k]
}

fc_log_p0_level <- ts(fc_log_p0_full[5:(4 + h)],
                      start = c(2024, 1), frequency = 4)
fc_log_p1_level <- ts(fc_log_p1_full[5:(4 + h)],
                      start = c(2024, 1), frequency = 4)

# True log levels in test period
test_log_p0 <- window(log_p0, start = c(2024, 1))
test_log_p1 <- window(log_p1, start = c(2024, 1))

cat("\n== TEST ACCURACY – Reconstructed LOG levels from VAR(Δ4 log + xreg) ==\n")
cat("\nlog p0:\n")
print(safe_accuracy(fc_log_p0_level, test_log_p0))
cat("\nlog p1:\n")
print(safe_accuracy(fc_log_p1_level, test_log_p1))

# Back-transform to LEVELS
fc_p0_level <- exp(fc_log_p0_level)
fc_p1_level <- exp(fc_log_p1_level)
fc_share_perc <- 100 * fc_p1_level / fc_p0_level

# True levels in test period
test_p0 <- window(ts_p0, start = c(2024, 1))
test_p1 <- window(ts_p1, start = c(2024, 1))
test_share <- window(ts_share, start = c(2024, 1))

cat("\n== TEST ACCURACY – Reconstructed LEVELS from VAR(Δ4 log + xreg) ==\n")
cat("\np0 level:\n")
print(safe_accuracy(fc_p0_level, test_p0))
cat("\np1 level:\n")
print(safe_accuracy(fc_p1_level, test_p1))
cat("\nshare (% of GDP):\n")
print(safe_accuracy(fc_share_perc, test_share))


# =========================
# Backtest Plots in Levels
# =========================
p_p0_level <- plot_backtest_line(
  ts_full   = ts_p0,
  yhat_test = fc_p0_level,
  start_yq  = c(2024, 1),
  title     = "VAR(Δ4 log + xreg) backtest – GDP p0 (level) 2024–2025",
  ylab      = "p0"
)

p_p1_level <- plot_backtest_line(
  ts_full   = ts_p1,
  yhat_test = fc_p1_level,
  start_yq  = c(2024, 1),
  title     = "VAR(Δ4 log + xreg) backtest – Agriculture p1 (level) 2024–2025",
  ylab      = "p1"
)

p_share_level <- plot_backtest_line(
  ts_full   = ts_share,
  yhat_test = fc_share_perc,
  start_yq  = c(2024, 1),
  title     = "VAR(Δ4 log + xreg) backtest – Agriculture share of GDP (%) 2024–2025",
  ylab      = "Share (%)"
)

# Print plots if running interactively
# p_p0_level; p_p1_level; p_share_level
