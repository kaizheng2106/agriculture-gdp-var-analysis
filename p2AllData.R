############################################################
# 0. Libraries & helpers
############################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(forecast)
  library(tseries)
  library(ggplot2)
  library(vars)
  library(urca)
  library(knitr)
})

# Accuracy helper (kept in case you need it later)
safe_accuracy <- function(fc, ytest) {
  A <- accuracy(fc, ytest)
  A[, intersect(c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1","Theil's U"),
                colnames(A)), drop = FALSE]
}

# ADF lag rule-of-thumb
adf_k <- function(x) trunc((length(na.omit(x)) - 1)^(1/3))


############################################################
# 1. Load & clean data
############################################################
df <- read.csv("original.csv", stringsAsFactors = FALSE)

norm_str <- function(x) tolower(trimws(as.character(x)))
df <- df %>%
  mutate(
    series = norm_str(series),
    sector = norm_str(sector),
    value  = as.numeric(gsub(",", "", as.character(value)))
  )

# Try several date formats
d0 <- suppressWarnings(dmy(df$date))
d1 <- suppressWarnings(ymd(df$date))
d2 <- suppressWarnings(mdy(df$date))
df$date <- coalesce(d0, d1, d2)
stopifnot(!all(is.na(df$date)))

abs_df <- df %>%
  filter(series == "abs", !is.na(value)) %>%
  arrange(date)

# GDP total (p0)
gdp  <- abs_df %>%
  dplyr::filter(sector %in% c("p0","gdp","gross domestic product","jumlah","kdnk")) %>%
  dplyr::arrange(date) %>%
  dplyr::select(date, p0 = value)

# Agriculture (p1)
agri <- abs_df %>%
  dplyr::filter(sector %in% c("p1","agriculture","agriculture, forestry & fishing","pertanian")) %>%
  dplyr::arrange(date) %>%
  dplyr::select(date, p1 = value)

dat <- inner_join(gdp, agri, by = "date") %>%
  arrange(date)


############################################################
# 2. Build quarterly ts objects (FULL SAMPLE)
############################################################
ts_p0 <- ts(dat$p0,
            start = c(year(min(dat$date)), quarter(min(dat$date))),
            frequency = 4)
ts_p1 <- ts(dat$p1,
            start = start(ts_p0),
            frequency = 4)

# Optional: agriculture share
ts_share <- ts(100 * dat$p1 / dat$p0,
               start = start(ts_p0),
               frequency = 4)

tq <- time(ts_p0)  # time index


############################################################
# 3. Deterministic terms: trend, seasonal dummies, pulses
############################################################

# Linear time trend (full sample)
trend_full <- 1:length(ts_p0)

# Quarterly seasonal dummies: 3 columns (Q1, Q2, Q3; Q4 baseline)
season_full <- seasonaldummy(ts_p0)

# Pulses: COVID and 2022Q2 shock
X_pulses <- cbind(
  covid_q2     = as.numeric(tq == 2020 + 1/4),
  covid_q3     = as.numeric(tq == 2020 + 2/4),
  shock_2022q2 = as.numeric(tq == 2022 + 1/4)
)

# Full exogenous matrix
exog_full <- cbind(trend = trend_full,
                   season_full,
                   X_pulses)

# Logs for cointegration & VAR
log_p0 <- log(ts_p0)
log_p1 <- log(ts_p1)

y_levels <- cbind(log_p0, log_p1)
colnames(y_levels) <- c("log_p0", "log_p1")


############################################################
# 4. Basic plots: levels vs diffs (for thesis)
############################################################

# Convert to data frames for ggplot
df_ts <- tibble(
  date = dat$date,
  log_p0 = as.numeric(log_p0),
  log_p1 = as.numeric(log_p1)
)

df_diff <- tibble(
  date = dat$date[-1],
  d1_log_p0 = diff(log_p0),
  d1_log_p1 = diff(log_p1)
)

df_sdiff <- tibble(
  date = dat$date[-(1:4)],
  d4_log_p0 = diff(log_p0, lag = 4),
  d4_log_p1 = diff(log_p1, lag = 4)
)

# Levels
ggplot(df_ts, aes(x = date)) +
  geom_line(aes(y = log_p0, colour = "log_p0")) +
  geom_line(aes(y = log_p1, colour = "log_p1")) +
  labs(title = "Log GDP (p0) and Log Agriculture GDP (p1)",
       y = "log value", colour = "") +
  theme_minimal()

# First differences
ggplot(df_diff, aes(x = date)) +
  geom_line(aes(y = d1_log_p0, colour = "Δ log_p0")) +
  geom_line(aes(y = d1_log_p1, colour = "Δ log_p1")) +
  labs(title = "First Differences of Logs",
       y = "Δ log", colour = "") +
  theme_minimal()

# Seasonal (lag-4) differences
ggplot(df_sdiff, aes(x = date)) +
  geom_line(aes(y = d4_log_p0, colour = "Δ4 log_p0")) +
  geom_line(aes(y = d4_log_p1, colour = "Δ4 log_p1")) +
  labs(title = "Seasonal (lag-4) Differences of Logs",
       y = "Δ4 log", colour = "") +
  theme_minimal()


############################################################
# 5. Unit-root tests: ADF + KPSS (FULL SAMPLE)
############################################################

# Differences (full sample)
d1_log_p0 <- diff(log_p0, lag = 1)
d1_log_p1 <- diff(log_p1, lag = 1)
d4_log_p0 <- diff(log_p0, lag = 4)
d4_log_p1 <- diff(log_p1, lag = 4)

# ADF on levels
adf_logp0_lvl <- adf.test(log_p0, k = adf_k(log_p0))
adf_logp1_lvl <- adf.test(log_p1, k = adf_k(log_p1))

# ADF on first diffs
adf_logp0_diff <- adf.test(d1_log_p0, k = adf_k(d1_log_p0))
adf_logp1_diff <- adf.test(d1_log_p1, k = adf_k(d1_log_p1))

# ADF on seasonal diffs
adf_logp0_sdiff <- adf.test(d4_log_p0, k = adf_k(d4_log_p0))
adf_logp1_sdiff <- adf.test(d4_log_p1, k = adf_k(d4_log_p1))

# KPSS (level-stationary) on diffs
kpss_d1_p0 <- kpss.test(d1_log_p0, null = "Level")
kpss_d1_p1 <- kpss.test(d1_log_p1, null = "Level")
kpss_d4_p0 <- kpss.test(d4_log_p0, null = "Level")
kpss_d4_p1 <- kpss.test(d4_log_p1, null = "Level")

# Print summaries
adf_logp0_lvl
adf_logp1_lvl
adf_logp0_diff
adf_logp1_diff
adf_logp0_sdiff
adf_logp1_sdiff

kpss_d1_p0
kpss_d1_p1
kpss_d4_p0
kpss_d4_p1

# Summary table for thesis (FULL SAMPLE)
adf_kpss_table <- tibble(
  Series = c("log_p0", "log_p1",
             "Δ log_p0", "Δ log_p1",
             "Δ4 log_p0", "Δ4 log_p1"),
  ADF_stat = c(adf_logp0_lvl$statistic,
               adf_logp1_lvl$statistic,
               adf_logp0_diff$statistic,
               adf_logp1_diff$statistic,
               adf_logp0_sdiff$statistic,
               adf_logp1_sdiff$statistic),
  ADF_p    = c(adf_logp0_lvl$p.value,
               adf_logp1_lvl$p.value,
               adf_logp0_diff$p.value,
               adf_logp1_diff$p.value,
               adf_logp0_sdiff$p.value,
               adf_logp1_sdiff$p.value),
  KPSS_stat = c(NA, NA,
                kpss_d1_p0$statistic,
                kpss_d1_p1$statistic,
                kpss_d4_p0$statistic,
                kpss_d4_p1$statistic),
  KPSS_p    = c(NA, NA,
                kpss_d1_p0$p.value,
                kpss_d1_p1$p.value,
                kpss_d4_p0$p.value,
                kpss_d4_p1$p.value)
)

kable(adf_kpss_table, digits = 4,
      caption = "ADF and KPSS tests on levels and differences (full sample)")


############################################################
# 6. Trend vs stochastic trend: residual ADF & KPSS on levels
############################################################

t_full <- 1:length(log_p0)
lm_p0 <- lm(log_p0 ~ t_full)
lm_p1 <- lm(log_p1 ~ t_full)

summary(lm_p0)
summary(lm_p1)

res_p0 <- resid(lm_p0)
res_p1 <- resid(lm_p1)

# ADF on residuals (test for trend-stationarity)
adf_res_p0 <- adf.test(res_p0, k = adf_k(res_p0))
adf_res_p1 <- adf.test(res_p1, k = adf_k(res_p1))
adf_res_p0
adf_res_p1

# KPSS for trend-stationarity (on levels)
kpss_logp0_mu    <- kpss.test(log_p0, null = "Level")
kpss_logp0_trend <- kpss.test(log_p0, null = "Trend")
kpss_logp1_mu    <- kpss.test(log_p1, null = "Level")
kpss_logp1_trend <- kpss.test(log_p1, null = "Trend")

kpss_logp0_mu
kpss_logp0_trend
kpss_logp1_mu
kpss_logp1_trend

# Plot residuals (for thesis)
res_df <- tibble(
  date = dat$date,
  res_p0 = res_p0,
  res_p1 = res_p1
)

ggplot(res_df, aes(x = date)) +
  geom_line(aes(y = res_p0, colour = "res(log_p0)")) +
  geom_line(aes(y = res_p1, colour = "res(log_p1)")) +
  labs(title = "Residuals after removing linear trend",
       y = "residual", colour = "") +
  theme_minimal()


############################################################
# 7. Johansen cointegration test with seasonal dummies
############################################################

S_full <- season_full  # seasonal dummies for full sample

# Choose lag order in levels (including seasonal dummies as exogenous)
lag_sel <- VARselect(y_levels,
                     lag.max = 4,
                     type    = "const",
                     exogen  = S_full)
lag_sel$selection
p_opt <- lag_sel$selection["AIC(n)"]
K_joh <- max(2, as.integer(p_opt) + 1)  # Johansen uses K-1 lags in Δ

# Johansen with trend in cointegration + seasonal dummies
joh_trend <- ca.jo(
  y_levels,
  type   = "trace",
  ecdet  = "trend",    # trend in cointegration space
  K      = K_joh,
  dumvar = S_full      # seasonal dummies only
)
summary(joh_trend)

# Engle-Granger-style residual test (validation)
trend_ts <- ts(1:nrow(y_levels), start = start(y_levels), frequency = 4)
Y_aug <- cbind(
  log_p0 = y_levels[,1],
  log_p1 = y_levels[,2],
  trend  = trend_ts
)

beta_vec <- joh_trend@V[,1]   # first cointegration vector
ect      <- Y_aug %*% beta_vec
ect      <- ts(ect, start = start(y_levels), frequency = 4)

autoplot(ect) +
  ggtitle("Cointegration Residual (ECT) from Johansen") +
  ylab("ECT") +
  theme_minimal()

# ADF on ECT: MUST be stationary if cointegrated
adf_ect <- adf.test(as.numeric(ect), k = adf_k(ect))
adf_ect

# If ADF on ect fails to reject unit root => no valid cointegration => use VAR in diffs


############################################################
# 8. FINAL MODELS: VAR on stationary transforms (FULL SAMPLE)
############################################################
# Model A: VAR on first differences with seasonal dummies + pulses
# Model B: VAR on seasonal differences with pulses only

############################################################
# 8A. Model A: VAR in Δ log with seasonal dummies + pulses
############################################################

d1_log_p0 <- diff(log_p0)
d1_log_p1 <- diff(log_p1)
y_d1 <- cbind(d1_log_p0, d1_log_p1)
colnames(y_d1) <- c("d1_log_p0", "d1_log_p1")

# Align exogenous for differenced series (drop first obs)
S_d1      <- season_full[-1, , drop = FALSE]
X_pulses1 <- X_pulses[-1, , drop = FALSE]
X_d1      <- cbind(S_d1, X_pulses1)

# Lag selection
lag_sel_d1 <- VARselect(y_d1, lag.max = 4, type = "const", exogen = X_d1)
lag_sel_d1$selection
p_d1 <- as.integer(lag_sel_d1$selection["AIC(n)"])

var_d1 <- VAR(
  y_d1,
  p      = p_d1,
  type   = "const",
  exogen = X_d1
)
summary(var_d1)

# Residual diagnostics
serial.test(var_d1, lags.pt = 8, type = "PT.asymptotic")  # Portmanteau
arch.test(var_d1, lags.multi = 4)
normality.test(var_d1)

# Optional: plot residuals
plot(resid(var_d1), main = "VAR(Δ log) residuals")


############################################################
# 8B. Model B: VAR in Δ4 log with pulses only
############################################################

d4_log_p0 <- diff(log_p0, lag = 4)
d4_log_p1 <- diff(log_p1, lag = 4)
y_d4 <- cbind(d4_log_p0, d4_log_p1)
colnames(y_d4) <- c("d4_log_p0", "d4_log_p1")

# Align pulses with Δ4 series using time index
t_y_d4 <- time(y_d4)
X_p_d4 <- X_pulses[match(t_y_d4, tq), , drop = FALSE]

# Lag selection
lag_sel_d4 <- VARselect(y_d4, lag.max = 4, type = "const", exogen = X_p_d4)
lag_sel_d4$selection
p_d4 <- as.integer(lag_sel_d4$selection["AIC(n)"])

var_d4 <- VAR(
  y_d4,
  p      = p_d4,
  type   = "const",
  exogen = X_p_d4
)
summary(var_d4)

# Residual diagnostics
serial.test(var_d4, lags.pt = 8, type = "PT.asymptotic")
arch.test(var_d4, lags.multi = 4)
normality.test(var_d4)

plot(resid(var_d4), main = "VAR(Δ4 log) residuals")


############################################################
# 9. Residual ACFs and BG tests (FULL SAMPLE)
############################################################

res_d1 <- resid(var_d1)
par(mfrow = c(2, 1))
acf(res_d1[, "d1_log_p0"], main = "ACF residuals Δlog_p0 (Model A)")
acf(res_d1[, "d1_log_p1"], main = "ACF residuals Δlog_p1 (Model A)")

res_d4 <- resid(var_d4)
par(mfrow = c(2, 1))
acf(res_d4[, "d4_log_p0"], main = "ACF residuals Δ4log_p0 (Model B)")
acf(res_d4[, "d4_log_p1"], main = "ACF residuals Δ4log_p1 (Model B)")
par(mfrow = c(1, 1))

serial.test(var_d1, type = "BG", lags.bg = 4)
serial.test(var_d4, type = "BG", lags.bg = 4)


############################################################
# END OF SCRIPT (FULL SAMPLE VERSION)
############################################################

library(vars)
roots(var_d1)   # var_model is the output from VAR()
roots(var_d4)   # var_model is the output from VAR()
