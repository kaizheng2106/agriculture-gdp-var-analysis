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
  library(uroot)
  library(knitr)
  library(quantmod)
})

# Accuracy helper (same as you had)
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

agri <- abs_df %>%
  dplyr::filter(sector %in% c("p1","agriculture","agriculture, forestry & fishing","pertanian")) %>%
  dplyr::arrange(date) %>%
  dplyr::select(date, p1 = value)

dat <- inner_join(gdp, agri, by = "date") %>%
  arrange(date)

############################################################
# 2. Build quarterly ts objects
############################################################
ts_p0 <- ts(dat$p0,
            start = c(year(min(dat$date)), quarter(min(dat$date))),
            frequency = 4)
ts_p1 <- ts(dat$p1,
            start = start(ts_p0),
            frequency = 4)

log_p0 <- log(ts_p0)
log_p1 <- log(ts_p1)

# 5. Perform STL Decomposition
# s.window = "periodic" forces the seasonal pattern to be constant over time

stl_p0 <- stl(ts_p0, s.window = "periodic")
stl_p1 <- stl(ts_p1, s.window = "periodic")
stl_lp0 <- stl(log_p0, s.window = "periodic")
stl_lp1 <- stl(log_p1, s.window = "periodic")

# 6. Plot the Results
# Plot for log_p0
plot(stl_p0, main = "STL Decomposition of p0")

# Plot for log_p1
plot(stl_p1, main = "STL Decomposition of p1")

# Plot for log_p0
plot(stl_lp0, main = "STL Decomposition of log_p0")

# Plot for log_p1
plot(stl_lp1, main = "STL Decomposition of log_p1")

# Optional: agriculture share
ts_share <- ts(100 * dat$p1 / dat$p0,
               start = start(ts_p0),
               frequency = 4)

tq <- time(ts_p0)  # time index


############################################################
# 3. Deterministic terms: trend, seasonal dummies, pulses
############################################################

# Linear time trend (for full sample)
trend_full <- 1:length(ts_p0)

# Quarterly seasonal dummies: 3 columns (Q1, Q2, Q3; Q4 is baseline)
season_full <- seasonaldummy(ts_p0)

# Pulses: COVID and 2022Q2 shock
X_pulses <- cbind(
  covid_q2     = as.numeric(tq == 2020 + 1/4),
  covid_q3     = as.numeric(tq == 2020 + 2/4),
  shock_2022q2 = as.numeric(tq == 2022 + 1/4)
)

# For convenience, full exogenous matrix
exog_full <- cbind(trend = trend_full,
                   season_full,
                   X_pulses)

# Logs for cointegration
log_p0 <- log(ts_p0)
log_p1 <- log(ts_p1)

y_levels <- cbind(log_p0, log_p1)
colnames(y_levels) <- c("log_p0", "log_p1")


############################################################
# 4. Train / test split (up to 2023Q4 as train)
############################################################
cut_year <- c(2023, 4)
cut_time <- 2023 + 3/4  # 2023Q4 in decimal

y_train <- window(y_levels, end   = cut_year)
y_test  <- window(y_levels, start = c(2024, 1))

train_mask <- (time(ts_p0) <= cut_time)
test_mask  <- !train_mask

exog_train <- exog_full[train_mask, , drop = FALSE]
exog_test  <- exog_full[test_mask,  , drop = FALSE]

stopifnot(nrow(y_train) == nrow(exog_train))


############################################################
# 5. Basic plots: levels vs diffs (for thesis)
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
# 6. Unit-root tests: ADF + KPSS (TRAIN sample)
############################################################

# TRAIN sample series
ts_p0_tr <- window(ts_p0, end = cut_year)
ts_p1_tr <- window(ts_p1, end = cut_year)
log_p0_tr <- window(log_p0, end = cut_year)
log_p1_tr <- window(log_p1, end = cut_year)

# Differences (train)
d1_log_p0_tr <- diff(log_p0_tr, lag = 1)
d1_log_p1_tr <- diff(log_p1_tr, lag = 1)
d4_log_p0_tr <- diff(log_p0_tr, lag = 4)
d4_log_p1_tr <- diff(log_p1_tr, lag = 4)

# ADF on levels
adf_logp0_lvl_tr <- adf.test(log_p0_tr, k = adf_k(log_p0_tr))
adf_logp1_lvl_tr <- adf.test(log_p1_tr, k = adf_k(log_p1_tr))

# ADF on first diffs
adf_logp0_diff_tr <- adf.test(d1_log_p0_tr, k = adf_k(d1_log_p0_tr))
adf_logp1_diff_tr <- adf.test(d1_log_p1_tr, k = adf_k(d1_log_p1_tr))

# ADF on seasonal diffs
adf_logp0_sdiff_tr <- adf.test(d4_log_p0_tr, k = adf_k(d4_log_p0_tr))
adf_logp1_sdiff_tr <- adf.test(d4_log_p1_tr, k = adf_k(d4_log_p1_tr))

# KPSS (level-stationary) on diffs
kpss_d1_p0 <- kpss.test(d1_log_p0_tr, null = "Level")
kpss_d1_p1 <- kpss.test(d1_log_p1_tr, null = "Level")
kpss_d4_p0 <- kpss.test(d4_log_p0_tr, null = "Level")
kpss_d4_p1 <- kpss.test(d4_log_p1_tr, null = "Level")

# Print summaries
adf_logp0_lvl_tr
adf_logp1_lvl_tr
adf_logp0_diff_tr
adf_logp1_diff_tr
adf_logp0_sdiff_tr
adf_logp1_sdiff_tr

kpss_d1_p0
kpss_d1_p1
kpss_d4_p0
kpss_d4_p1

# (Optional) small summary table for thesis
adf_kpss_table <- tibble(
  Series = c("log_p0_tr", "log_p1_tr",
             "Δ log_p0_tr", "Δ log_p1_tr",
             "Δ4 log_p0_tr", "Δ4 log_p1_tr"),
  ADF_stat = c(adf_logp0_lvl_tr$statistic,
               adf_logp1_lvl_tr$statistic,
               adf_logp0_diff_tr$statistic,
               adf_logp1_diff_tr$statistic,
               adf_logp0_sdiff_tr$statistic,
               adf_logp1_sdiff_tr$statistic),
  ADF_p    = c(adf_logp0_lvl_tr$p.value,
               adf_logp1_lvl_tr$p.value,
               adf_logp0_diff_tr$p.value,
               adf_logp1_diff_tr$p.value,
               adf_logp0_sdiff_tr$p.value,
               adf_logp1_sdiff_tr$p.value),
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
      caption = "ADF and KPSS tests on levels and differences (training sample)")


############################################################
# 7. Trend vs stochastic trend: residual ADF & KPSS on levels
############################################################

t_tr <- 1:length(log_p0_tr)
lm_p0 <- lm(log_p0_tr ~ t_tr)
lm_p1 <- lm(log_p1_tr ~ t_tr)

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
kpss_logp0_tr_mu   <- kpss.test(log_p0_tr, null = "Level")
kpss_logp0_tr_trend<- kpss.test(log_p0_tr, null = "Trend")
kpss_logp1_tr_mu   <- kpss.test(log_p1_tr, null = "Level")
kpss_logp1_tr_trend<- kpss.test(log_p1_tr, null = "Trend")

kpss_logp0_tr_mu
kpss_logp0_tr_trend
kpss_logp1_tr_mu
kpss_logp1_tr_trend

# Plot residuals (for thesis)
res_df <- tibble(
  date = dat$date[train_mask],
  res_p0 = res_p0,
  res_p1 = res_p1
)

ggplot(res_df, aes(x = date)) +
  geom_line(aes(y = res_p0, colour = "res(log_p0_tr)")) +
  geom_line(aes(y = res_p1, colour = "res(log_p1_tr)")) +
  labs(title = "Residuals after removing linear trend",
       y = "residual", colour = "") +
  theme_minimal()

############################################################
# 7B. Trend & Seasonal vs Stochastic: Residual Analysis
############################################################

# 1. Prepare Regressors
# Time Trend
t_tr <- 1:length(log_p0_tr)
# Seasonal Dummies (for training set)
S_train <- window(seas, end = c(2023, 4)) 
# Convert to matrix for lm
S_mat <- as.matrix(S_train)

# 2. Run Regression: Log GDP ~ Trend + Seasonality
# This attempts to capture all deterministic movement
lm_p0_full <- lm(lp0_train ~ t_tr + S_mat)
lm_p1_full <- lm(lp1_train ~ t_tr + S_mat)

# 3. Extract Residuals ("The Noise")
res_p0_clean <- resid(lm_p0_full)
res_p1_clean <- resid(lm_p1_full)

# 4. ADF Test on these Residuals
# Note: We use type="none" because mean/trend are already removed.
# Critical Value Warning: The standard CV (-1.95) is technically too low 
# because we estimated parameters, but for visual proof, it's fine.
adf_res_p0 <- ur.df(res_p0_clean, type="none", lags=4, selectlags="AIC")
adf_res_p1 <- ur.df(res_p1_clean, type="none", lags=4, selectlags="AIC")

print("--- ADF on Detrended & Deseasonalized Residuals ---")
cat("P0 t-stat:", adf_res_p0@teststat[1], "\n")
cat("P1 t-stat:", adf_res_p1@teststat[1], "\n")

# 5. The Visual Proof (Most Important for Thesis)
# If these lines wander away from 0, you need Differencing.
res_df <- data.frame(
  date = as.numeric(time(lp0_train)),
  res_p0 = as.numeric(res_p0_clean),
  res_p1 = as.numeric(res_p1_clean)
)

ggplot(res_df, aes(x = date)) +
  geom_line(aes(y = res_p0, colour = "GDP (p0)")) +
  geom_line(aes(y = res_p1, colour = "Agri (p1)")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Residuals after removing Trend AND Seasonality",
       subtitle = "If these lines wander (are not flat noise), the series is I(1).",
       y = "Deviations from Trend", colour = "Series") +
  theme_minimal()

############################################################
# 8. Johansen cointegration test with seasonal dummies
############################################################

# Seasonal dummies for TRAIN only
S_train <- season_full[train_mask, , drop = FALSE]

# Choose lag order in levels (including seasonal dummies as exogenous)
lag_sel <- VARselect(y_train,
                     lag.max = 4,
                     type    = "trend",
                     exogen  = S_train)
lag_sel$selection
p_opt <- lag_sel$selection["AIC(n)"]
K_joh <- max(2, as.integer(p_opt))  # Johansen uses K-1 lags in Δ 

# Johansen with trend in cointegration + seasonal dummies
joh_trend <- ca.jo(
  y_train,
  type   = "trace",
  ecdet  = "trend",    # trend in cointegration space
  K      = K_joh,
  dumvar = S_train     # seasonal dummies only
)
summary(joh_trend)

# Engle-Granger style residual test (validation)
trend_ts <- ts(1:nrow(y_train), start = start(y_train), frequency = 4)
Y_aug <- cbind(
  log_p0 = y_train[,1],
  log_p1 = y_train[,2],
  trend  = trend_ts
)

beta_vec <- joh_trend@V[,1]   # first cointegration vector
ect      <- Y_aug %*% beta_vec
ect      <- ts(ect, start = start(y_train), frequency = 4)

autoplot(ect) +
  ggtitle("Cointegration Residual (ECT) from Johansen") +
  ylab("ECT") +
  theme_minimal()

# ADF on ECT: MUST be stationary if cointegrated
adf_ect <- adf.test(as.numeric(ect), k = adf_k(ect))
adf_ect

# If ADF on ect fails to reject unit root => no valid cointegration => use VAR in diffs

############################################################
# 8B. Johansen Cointegration Test (Robust with Pulses)
############################################################
library(vars)
library(urca)

# --- 1. Prepare Exogenous Variables (Train Only) ---
# We need to combine Seasonal Dummies + Pulses for the training period

# Subset Seasonals
S_train <- window(season_full, end = cut_year)

# Subset Pulses
X_pulses_train <- window(X_pulses, end = cut_year)

# Combine them into one matrix for the test
# dumvar requires a matrix, not a ts object sometimes, so we convert.
exo_train <- cbind(as.matrix(S_train), as.matrix(X_pulses_train))
colnames(exo_train) <- c("S1","S2","S3", "CovidQ2", "CovidQ3", "Shock22")

# --- 2. Lag Selection ---
# We select lags for the VAR in Levels
lag_sel <- VARselect(y_train, 
                     lag.max = 4, 
                     type = "trend", 
                     exogen = exo_train) # Important: Include pulses here!

p_opt <- lag_sel$selection["AIC(n)"]
cat("Optimal Lag (AIC):", p_opt, "\n")

# Johansen uses K = p_opt (lags in levels). 
# Note: Sometimes defined as K-1 in diffs. ca.jo takes the 'levels' lag K.
K_joh <- max(2, as.integer(p_opt))

# --- 3. Run Johansen Test (Trace) ---
# We use ecdet="const" (Intercept in cointegration). 
# "trend" assumes the GAP between GDP and Agri grows linearly forever, 
# which is a very strong assumption. "const" is safer.

joh_test <- ca.jo(
  y_train,
  type = "trace",
  ecdet = "const",      # Try "const" first. 
  K = K_joh,
  dumvar = exo_train    # Seasonality + Pulses included
)

summary(joh_test)

# --- 4. Validation: Check if the Residuals (ECT) are Stationary ---
# If variables are cointegrated, the linear combination (residuals) MUST be stationary.

# Extract the cointegration residuals (ECT) mathematically
# We use cajorls() to get the restricted VECM parameters
vecm_fit <- cajorls(joh_test, r = 1) # Assuming r=1 for visualization
ect_series <- vecm_fit$rlm$residuals[,1] # Take residuals from first equation

# Plot ECT
ts.plot(ect_series, main="Cointegration Residuals (ECT)", ylab="Value")
abline(h=0, col="red", lty=2)

# Run ADF on the ECT
# If this rejects Null (is Stationary), Cointegration is valid.
adf_ect <- ur.df(ect_series, type="none", lags=1, selectlags="AIC")
cat("ADF on ECT t-stat:", adf_ect@teststat[1], "\n")
cat("ADF 5% Critical Value:", adf_ect@cval[1,2], "\n")

############################################################
# 9. FINAL MODELS: VAR on stationary transforms
############################################################
# You can choose either:
#   Model A: VAR on first differences with seasonal dummies + pulses
#   Model B: VAR on seasonal differences with pulses only
# Both are shown; pick one as your main model.

############################################################
# 9A. Model A: VAR in Δ log with seasonal dummies + pulses
############################################################

d1_log_p0_tr <- diff(log_p0_tr)
d1_log_p1_tr <- diff(log_p1_tr)
y_d1 <- cbind(d1_log_p0_tr, d1_log_p1_tr)
colnames(y_d1) <- c("d1_log_p0", "d1_log_p1")

# Align exogenous for differenced series (drop first obs)
S_tr_d1 <- S_train[-1, , drop = FALSE]
X_pulses_tr <- X_pulses[train_mask, , drop = FALSE]
X_p_tr_d1 <- X_pulses_tr[-1, , drop = FALSE]
X_d1 <- cbind(S_tr_d1, X_p_tr_d1)

# Lag selection
lag_sel_d1 <- VARselect(y_d1, lag.max = 4, type = "trend", exogen = X_d1)
lag_sel_d1$selection
p_d1 <- as.integer(lag_sel_d1$selection["AIC(n)"])

var_d1 <- VAR(
  y_d1,
  p      = p_d1,
  type   = "trend",
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
# 9B. Model B: VAR in Δ4 log with pulses only
############################################################

d4_log_p0_tr <- diff(log_p0_tr, lag = 4)
d4_log_p1_tr <- diff(log_p1_tr, lag = 4)
y_d4 <- cbind(d4_log_p0_tr, d4_log_p1_tr)
colnames(y_d4) <- c("d4_log_p0", "d4_log_p1")

# Align pulses with Δ4 series (drop first 4 obs)
X_p_tr <- X_pulses[train_mask, , drop = FALSE]
t_y_d4 <- time(y_d4)
X_p_d4 <- X_p_tr[match(t_y_d4, tq[train_mask]), , drop = FALSE]

# Lag selection
lag_sel_d4 <- VARselect(y_d4, lag.max = 4, type = "trend", exogen = X_p_d4)
lag_sel_d4$selection
p_d4 <- as.integer(lag_sel_d4$selection["AIC(n)"])

var_d4 <- VAR(
  y_d4,
  p      = p_d4,
  type   = "trend",
  exogen = X_p_d4
)
summary(var_d4)

# Residual diagnostics
serial.test(var_d4, lags.pt = 8, type = "PT.asymptotic")
arch.test(var_d4, lags.multi = 4)
normality.test(var_d4)

plot(resid(var_d4), main = "VAR(Δ4 log) residuals")

############################################################
# 9C. Model C: VAR in Combined Diff (ΔΔ4 log)
############################################################

# 1. Create Double Differenced Series (Train)
# We take the seasonally differenced series (d4) and diff them again (d1)
d1d4_log_p0_tr <- diff(d4_log_p0_tr, lag = 1)
d1d4_log_p1_tr <- diff(d4_log_p1_tr, lag = 1)

y_d1d4 <- cbind(d1d4_log_p0_tr, d1d4_log_p1_tr)
colnames(y_d1d4) <- c("d1d4_log_p0", "d1d4_log_p1")

# 2. Align Pulses
# We lose 5 observations total (4 from d4, 1 from d1).
# d4_log_pX_tr already dropped 4. We just need to drop 1 more from that timeline.
# Matching time indices is the safest way.

t_y_d1d4 <- time(y_d1d4)
X_p_tr_all <- X_pulses[train_mask, , drop = FALSE]
X_p_d1d4   <- X_p_tr_all[match(t_y_d1d4, tq[train_mask]), , drop = FALSE]

# 3. Lag Selection
# Note: For double differenced data, usually 'const' (intercept) is zero, 
# but 'trend' is definitely gone. We use type="const" or "none". 
# Let's use "const" to be safe.
lag_sel_d1d4 <- VARselect(y_d1d4, lag.max = 4, type = "const", exogen = X_p_d1d4)
p_d1d4 <- as.integer(lag_sel_d1d4$selection["AIC(n)"])
cat("Optimal Lag for Model C:", p_d1d4, "\n")

# 4. Estimation
var_d1d4 <- VAR(
  y_d1d4,
  p      = p_d1d4,
  type   = "const", 
  exogen = X_p_d1d4
)
summary(var_d1d4)

# 5. Diagnostics
# 1. Serial Correlation (Portmanteau)
# Checks if there are patterns left in the noise
serial.test(var_d1d4, lags.pt = 8, type = "PT.asymptotic")

# 2. Heteroscedasticity (ARCH)
# Checks if the volatility (variance) changes over time
arch.test(var_d1d4, lags.multi = 4)

# 3. Normality (Jarque-Bera)
# Checks if the residuals follow a Bell Curve
norm_test_c <- normality.test(var_d1d4)
print(norm_test_c$JB)       # Overall test
print(norm_test_c$Skewness) # Symmetry
print(norm_test_c$Kurtosis) # Fat tails

# Stability (Roots)
roots_c <- roots(var_d1d4)
cat("Roots for Model C (modulus < 1 is stable):\n")
print(roots_c)

# Plot Residuals
plot(resid(var_d1d4), main = "VAR(ΔΔ4 log) residuals (Model C)")
acf(resid(var_d1d4))

############################################################
# 10. Out-of-sample forecasting (train-only estimation)
############################################################

# 10.1 Build full-sample differenced series ----------------

# First differences of logs for the whole sample
d1_log_p0_all <- diff(log_p0)
d1_log_p1_all <- diff(log_p1)
y_d1_all <- cbind(d1_log_p0_all, d1_log_p1_all)
colnames(y_d1_all) <- c("d1_log_p0", "d1_log_p1")

# Seasonal (lag‑4) differences of logs for the whole sample
d4_log_p0_all <- diff(log_p0, lag = 4)
d4_log_p1_all <- diff(log_p1, lag = 4)
y_d4_all <- cbind(d4_log_p0_all, d4_log_p1_all)
colnames(y_d4_all) <- c("d4_log_p0", "d4_log_p1")

# 10.2 Extract TEST parts of differenced series ------------
# Test sample starts at 2024Q1 (consistent with y_test above)
y_d1_test <- window(y_d1_all, start = c(2024, 1))  # Δlog test
y_d4_test <- window(y_d4_all, start = c(2024, 1))  # Δ4log test

h_d1 <- if (length(y_d1_test) == 0) 0 else nrow(as.matrix(y_d1_test))
h_d4 <- if (length(y_d4_test) == 0) 0 else nrow(as.matrix(y_d4_test))

cat("Horizon: Model A (Δlog) =", h_d1,
    "quarters; Model B (Δ4log) =", h_d4, "quarters\n")

# 10.3 Exogenous matrices for forecast horizon -------------

# Seasonal dummies and pulses for TEST in levels
S_test        <- season_full[test_mask, , drop = FALSE]
X_pulses_test <- X_pulses[test_mask, , drop = FALSE]

# For Model A (Δlog), estimation used X_d1 = cbind(S_tr_d1, X_p_tr_d1).
# For forecasting, use the same columns for the test horizon.
if (h_d1 > 0) {
  X_d1_fore <- cbind(S_test, X_pulses_test)
  
  # keep only the columns that were used in estimation, in the same order
  X_d1_fore <- X_d1_fore[seq_len(h_d1), colnames(X_d1), drop = FALSE]
}

# For Model B (Δ4log), only the pulse dummies are used as exogenous (X_p_d4).
# Build a pulse matrix for the y_d4_test time index and match colnames(X_p_d4).
if (h_d4 > 0) {
  t_test_d4 <- time(y_d4_test)
  
  X_p_d4_test <- cbind(
    covid_q2     = as.numeric(t_test_d4 == (2020 + 1/4)),
    covid_q3     = as.numeric(t_test_d4 == (2020 + 2/4)),
    shock_2022q2 = as.numeric(t_test_d4 == (2022 + 1/4))
  )
  
  # Ensure same column ordering as in estimation
  X_p_d4_test <- X_p_d4_test[, colnames(X_p_d4), drop = FALSE]
}

# 10.4 Model A: VAR in Δlog --------------------------------

if (h_d1 > 0) {
  
  predA <- predict(var_d1, n.ahead = h_d1, dumvar = X_d1_fore)
  
  fc_d1_p0 <- as.numeric(predA$fcst[["d1_log_p0"]][, "fcst"])
  fc_d1_p1 <- as.numeric(predA$fcst[["d1_log_p1"]][, "fcst"])
  
  fc_d1_p0_ts <- ts(fc_d1_p0, start = start(y_d1_test), frequency = 4)
  fc_d1_p1_ts <- ts(fc_d1_p1, start = start(y_d1_test), frequency = 4)
  
  par(mfrow = c(2,1))
  ts.plot(y_d1_test[, "d1_log_p0"], fc_d1_p0_ts,
          col = c("black","red"), lty = c(1,2),
          main = "Model A: Δlog(p0) actual vs forecast",
          ylab = "Δlog p0")
  legend("topleft", c("actual","forecast"),
         col = c("black","red"), lty = c(1,2), bty = "n")
  
  ts.plot(y_d1_test[, "d1_log_p1"], fc_d1_p1_ts,
          col = c("black","blue"), lty = c(1,2),
          main = "Model A: Δlog(p1) actual vs forecast",
          ylab = "Δlog p1")
  legend("topleft", c("actual","forecast"),
         col = c("black","blue"), lty = c(1,2), bty = "n")
  par(mfrow = c(1,1))
  
  cat("\nModel A (Δlog) accuracy:\n")
  print(safe_accuracy(fc_d1_p0_ts, y_d1_test[, "d1_log_p0"]))
  print(safe_accuracy(fc_d1_p1_ts, y_d1_test[, "d1_log_p1"]))
  
  # Reconstruct log-level forecasts by cumulating Δlog
  last_log_p0 <- tail(window(log_p0, end = cut_year), 1)
  last_log_p1 <- tail(window(log_p1, end = cut_year), 1)
  
  log_p0_fc_A <- ts(cumsum(c(last_log_p0, fc_d1_p0))[-1],
                    start = start(y_d1_test), frequency = 4)
  log_p1_fc_A <- ts(cumsum(c(last_log_p1, fc_d1_p1))[-1],
                    start = start(y_d1_test), frequency = 4)
  
  if (!is.null(y_test) && length(y_test) > 0) {
    par(mfrow = c(2,1))
    ts.plot(y_test[, "log_p0"], log_p0_fc_A,
            col = c("black","red"), lty = c(1,2),
            main = "Model A: log p0 actual vs forecast",
            ylab = "log p0")
    legend("topleft", c("actual","forecast"),
           col = c("black","red"), lty = c(1,2), bty = "n")
    
    ts.plot(y_test[, "log_p1"], log_p1_fc_A,
            col = c("black","blue"), lty = c(1,2),
            main = "Model A: log p1 actual vs forecast",
            ylab = "log p1")
    legend("topleft", c("actual","forecast"),
           col = c("black","blue"), lty = c(1,2), bty = "n")
    par(mfrow = c(1,1))
  }
  
} else {
  cat("No Δlog test observations; skipping Model A.\n")
}

# 10.5 Model B: VAR in Δ4log -------------------------------

if (h_d4 > 0) {
  
  predB <- predict(var_d4, n.ahead = h_d4, dumvar = X_p_d4_test)
  
  fc_d4_p0 <- as.numeric(predB$fcst[["d4_log_p0"]][, "fcst"])
  fc_d4_p1 <- as.numeric(predB$fcst[["d4_log_p1"]][, "fcst"])
  
  fc_d4_p0_ts <- ts(fc_d4_p0, start = start(y_d4_test), frequency = 4)
  fc_d4_p1_ts <- ts(fc_d4_p1, start = start(y_d4_test), frequency = 4)
  
  par(mfrow = c(2,1))
  ts.plot(y_d4_test[, "d4_log_p0"], fc_d4_p0_ts,
          col = c("black","red"), lty = c(1,2),
          main = "Model B: Δ4log(p0) actual vs forecast",
          ylab = "Δ4log p0")
  legend("topleft", c("actual","forecast"),
         col = c("black","red"), lty = c(1,2), bty = "n")
  
  ts.plot(y_d4_test[, "d4_log_p1"], fc_d4_p1_ts,
          col = c("black","blue"), lty = c(1,2),
          main = "Model B: Δ4log(p1) actual vs forecast",
          ylab = "Δ4log p1")
  legend("topleft", c("actual","forecast"),
         col = c("black","blue"), lty = c(1,2), bty = "n")
  par(mfrow = c(1,1))
  
  cat("\nModel B (Δ4log) accuracy:\n")
  print(safe_accuracy(fc_d4_p0_ts, y_d4_test[, "d4_log_p0"]))
  print(safe_accuracy(fc_d4_p1_ts, y_d4_test[, "d4_log_p1"]))
  
  # Reconstruct log levels from Δ4log using last 4 in-sample logs
  last4_log_p0 <- as.numeric(tail(window(log_p0, end = cut_year), 4))
  last4_log_p1 <- as.numeric(tail(window(log_p1, end = cut_year), 4))
  
  buf_p0 <- c(last4_log_p0, rep(NA_real_, h_d4))
  buf_p1 <- c(last4_log_p1, rep(NA_real_, h_d4))
  
  for (i in seq_len(h_d4)) {
    buf_p0[4 + i] <- fc_d4_p0[i] + buf_p0[i]
    buf_p1[4 + i] <- fc_d4_p1[i] + buf_p1[i]
  }
  
  log_p0_fc_B <- ts(buf_p0[5:(4 + h_d4)],
                    start = start(y_d4_test), frequency = 4)
  log_p1_fc_B <- ts(buf_p1[5:(4 + h_d4)],
                    start = start(y_d4_test), frequency = 4)
  
  if (!is.null(y_test) && length(y_test) > 0) {
    par(mfrow = c(2,1))
    ts.plot(y_test[, "log_p0"], log_p0_fc_B,
            col = c("black","red"), lty = c(1,2),
            main = "Model B: log p0 actual vs forecast",
            ylab = "log p0")
    legend("topleft", c("actual","forecast"),
           col = c("black","red"), lty = c(1,2), bty = "n")
    
    ts.plot(y_test[, "log_p1"], log_p1_fc_B,
            col = c("black","blue"), lty = c(1,2),
            main = "Model B: log p1 actual vs forecast",
            ylab = "log p1")
    legend("topleft", c("actual","forecast"),
           col = c("black","blue"), lty = c(1,2), bty = "n")
    par(mfrow = c(1,1))
  }
  
} else {
  cat("No Δ4log test observations; skipping Model B.\n")
}

cat("\nSection 10 complete: forecasts produced from models estimated on training data only.\n")



#transform

## ---- Model A: transform log forecasts back to levels ----
# Actual GDP levels for the test period
p0_test_A <- window(ts_p0,
                    start = start(log_p0_fc_A),
                    end   = end(log_p0_fc_A))
p1_test_A <- window(ts_p1,
                    start = start(log_p1_fc_A),
                    end   = end(log_p1_fc_A))

# Forecasted levels (invert logs)
p0_fc_A <- ts(exp(log_p0_fc_A),
              start = start(log_p0_fc_A),
              frequency = 4)
p1_fc_A <- ts(exp(log_p1_fc_A),
              start = start(log_p1_fc_A),
              frequency = 4)

# Plot level actual vs forecast (Model A)
par(mfrow = c(2,1))
ts.plot(p0_test_A, p0_fc_A,
        col = c("black", "red"), lty = c(1, 2),
        main = "Model A: GDP (p0) – level actual vs forecast",
        ylab = "GDP level (p0)")
legend("topleft", c("actual", "forecast"),
       col = c("black", "red"), lty = c(1, 2), bty = "n")

ts.plot(p1_test_A, p1_fc_A,
        col = c("black", "blue"), lty = c(1, 2),
        main = "Model A: Agriculture GDP (p1) – level actual vs forecast",
        ylab = "GDP level (p1)")
legend("topleft", c("actual", "forecast"),
       col = c("black", "blue"), lty = c(1, 2), bty = "n")
par(mfrow = c(1,1))

# Optional: accuracy in levels
cat("\nModel A (LEVELS) accuracy:\n")
cat("p0:\n")
print(safe_accuracy(p0_fc_A, p0_test_A))
cat("p1:\n")
print(safe_accuracy(p1_fc_A, p1_test_A))


## ---- Model B: transform log forecasts back to levels ----
# Actual GDP levels for the test period
p0_test_B <- window(ts_p0,
                    start = start(log_p0_fc_B),
                    end   = end(log_p0_fc_B))
p1_test_B <- window(ts_p1,
                    start = start(log_p1_fc_B),
                    end   = end(log_p1_fc_B))

# Forecasted levels (invert logs)
p0_fc_B <- ts(exp(log_p0_fc_B),
              start = start(log_p0_fc_B),
              frequency = 4)
p1_fc_B <- ts(exp(log_p1_fc_B),
              start = start(log_p1_fc_B),
              frequency = 4)

# Plot level actual vs forecast (Model B)
par(mfrow = c(2,1))
ts.plot(p0_test_B, p0_fc_B,
        col = c("black", "red"), lty = c(1, 2),
        main = "Model B: GDP (p0) – level actual vs forecast",
        ylab = "GDP level (p0)")
legend("topleft", c("actual", "forecast"),
       col = c("black", "red"), lty = c(1, 2), bty = "n")

ts.plot(p1_test_B, p1_fc_B,
        col = c("black", "blue"), lty = c(1, 2),
        main = "Model B: Agriculture GDP (p1) – level actual vs forecast",
        ylab = "GDP level (p1)")
legend("topleft", c("actual", "forecast"),
       col = c("black", "blue"), lty = c(1, 2), bty = "n")
par(mfrow = c(1,1))

# Optional: accuracy in levels
cat("\nModel B (LEVELS) accuracy:\n")
cat("p0:\n")
print(safe_accuracy(p0_fc_B, p0_test_B))
cat("p1:\n")
print(safe_accuracy(p1_fc_B, p1_test_B))



# You can then reconstruct level forecasts if needed by summing back Δ4 logs,
# or just evaluate forecast accuracy on Δ4 logs if you transform the test data.

############################################################
# 10.6 Model C: Forecasting & Reconstruction
############################################################

# 1. Prepare Test Data for Evaluation
# Construct actual ΔΔ4 series for test set
d1d4_log_p0_all <- diff(d4_log_p0_all, lag=1)
d1d4_log_p1_all <- diff(d4_log_p1_all, lag=1)

y_d1d4_test <- window(cbind(d1d4_log_p0_all, d1d4_log_p1_all), start = c(2024, 1))
h_c <- nrow(y_d1d4_test)

if(h_c > 0) {
  # 2. Exogenous variables for Test
  # Pulses only (aligned to test time)
  t_test_c <- time(y_d1d4_test)
  X_p_test_c <- cbind(
    covid_q2     = as.numeric(t_test_c == (2020 + 1/4)),
    covid_q3     = as.numeric(t_test_c == (2020 + 2/4)),
    shock_2022q2 = as.numeric(t_test_c == (2022 + 1/4))
  )
  # Ensure column order matches estimation
  X_p_test_c <- X_p_test_c[, colnames(X_p_d1d4), drop = FALSE]
  
  # 3. Predict ΔΔ4
  predC <- predict(var_d1d4, n.ahead = h_c, dumvar = X_p_test_c)
  
  fc_d1d4_p0 <- as.numeric(predC$fcst[["d1d4_log_p0"]][, "fcst"])
  fc_d1d4_p1 <- as.numeric(predC$fcst[["d1d4_log_p1"]][, "fcst"])
  
  # 4. Reconstruction Step 1: Recover Δ4 (Seasonal Diffs)
  # Formula: Δ4_t = ΔΔ4_t + Δ4_{t-1}
  # We need the last observed Δ4 value from the Training set
  last_d4_p0 <- tail(d4_log_p0_tr, 1)
  last_d4_p1 <- tail(d4_log_p1_tr, 1)
  
  # Cumulate the changes
  fc_d4_p0_C <- ts(cumsum(c(last_d4_p0, fc_d1d4_p0))[-1], start=start(y_d1d4_test), frequency=4)
  fc_d4_p1_C <- ts(cumsum(c(last_d4_p1, fc_d1d4_p1))[-1], start=start(y_d1d4_test), frequency=4)
  
  # 5. Reconstruction Step 2: Recover Log Levels
  # Formula: Y_t = Δ4_t + Y_{t-4}
  # We need the last 4 observed Log Level values from the Training set
  last4_log_p0 <- as.numeric(tail(log_p0_tr, 4))
  last4_log_p1 <- as.numeric(tail(log_p1_tr, 4))
  
  buf_p0_c <- c(last4_log_p0, rep(NA_real_, h_c))
  buf_p1_c <- c(last4_log_p1, rep(NA_real_, h_c))
  
  # Iterative reconstruction
  for (i in seq_len(h_c)) {
    buf_p0_c[4 + i] <- fc_d4_p0_C[i] + buf_p0_c[i] # Current diff + value 4 quarters ago
    buf_p1_c[4 + i] <- fc_d4_p1_C[i] + buf_p1_c[i]
  }
  
  log_p0_fc_C <- ts(buf_p0_c[5:(4 + h_c)], start = start(y_d1d4_test), frequency = 4)
  log_p1_fc_C <- ts(buf_p1_c[5:(4 + h_c)], start = start(y_d1d4_test), frequency = 4)
  
  # 6. Plotting
  if (!is.null(y_test) && length(y_test) > 0) {
    par(mfrow = c(2,1))
    ts.plot(y_test[, "log_p0"], log_p0_fc_C,
            col = c("black","green"), lty = c(1,2),
            main = "Model C: log p0 actual vs forecast (ΔΔ4)",
            ylab = "log p0")
    legend("topleft", c("actual","forecast"), col = c("black","green"), lty = c(1,2), bty = "n")
    
    ts.plot(y_test[, "log_p1"], log_p1_fc_C,
            col = c("black","purple"), lty = c(1,2),
            main = "Model C: log p1 actual vs forecast (ΔΔ4)",
            ylab = "log p1")
    legend("topleft", c("actual","forecast"), col = c("black","purple"), lty = c(1,2), bty = "n")
    par(mfrow = c(1,1))
    
    # Levels Accuracy
    p0_fc_level_C <- exp(log_p0_fc_C)
    p1_fc_level_C <- exp(log_p1_fc_C)
    
    cat("\nModel C (LEVELS) accuracy:\n")
    print(safe_accuracy(p0_fc_level_C, exp(y_test[,"log_p0"])))
    print(safe_accuracy(p1_fc_level_C, exp(y_test[,"log_p1"])))
  }
  
} else {
  cat("No data for Model C forecast.\n")
}

###########################################################
# 11. IRF and FEVD – effect of agriculture (p1) on GDP (p0)
############################################################

# We work with the stationary VARs:
#   Model A: var_d1 (Δlog)
#   Model B: var_d4 (Δ4log)

############################################################
# 11A. IRF: shock to agriculture and response of GDP
############################################################

# Model A: impulse = Δlog_p1, response = Δlog_p0
irf_A_p1_to_p0 <- irf(
  var_d1,
  impulse  = "d1_log_p1",
  response = "d1_log_p0",
  n.ahead  = 12,      # horizons (quarters) – change if you want
  boot     = TRUE
)

# Model B: impulse = Δ4log_p1, response = Δ4log_p0
irf_B_p1_to_p0 <- irf(
  var_d4,
  impulse  = "d4_log_p1",
  response = "d4_log_p0",
  n.ahead  = 12,
  boot     = TRUE
)

# Plot IRFs
par(mfrow = c(1, 2))
plot(irf_A_p1_to_p0, main = "IRF (Model A): shock to Δlog(p1) on Δlog(p0)")
plot(irf_B_p1_to_p0, main = "IRF (Model B): shock to Δ4log(p1) on Δ4log(p0)")
par(mfrow = c(1, 1))


############################################################
# 11B. FEVD: share of p0 forecast variance explained by p1
############################################################

# Model A FEVD
fevd_A <- fevd(var_d1, n.ahead = 12)
# This is a list; the component for the p0 equation is:
#   fevd_A$d1_log_p0  (rows = horizons, cols = variables)

fevd_A_p0 <- fevd_A$d1_log_p0  # matrix: horizon × {d1_log_p0, d1_log_p1}

# Model B FEVD
fevd_B <- fevd(var_d4, n.ahead = 12)
fevd_B_p0 <- fevd_B$d4_log_p0  # horizon × {d4_log_p0, d4_log_p1}

# Simple plots of FEVD contributions of p1 to p0
par(mfrow = c(1, 2))
plot(1:12, fevd_A_p0[, "d1_log_p1"] * 100, type = "b",
     xlab = "Horizon (quarters)", ylab = "% variance",
     main = "FEVD (Model A): share of Δlog(p0) due to Δlog(p1)")
plot(1:12, fevd_B_p0[, "d4_log_p1"] * 100, type = "b",
     xlab = "Horizon (quarters)", ylab = "% variance",
     main = "FEVD (Model B): share of Δ4log(p0) due to Δ4log(p1)")
par(mfrow = c(1, 1))

# Optional: nice tables for thesis
fevd_table_A <- tibble(
  horizon = 1:12,
  `own shock (p0)` = fevd_A_p0[, "d1_log_p0"] * 100,
  `agri shock (p1)` = fevd_A_p0[, "d1_log_p1"] * 100
)

fevd_table_B <- tibble(
  horizon = 1:12,
  `own shock (p0)` = fevd_B_p0[, "d4_log_p0"] * 100,
  `agri shock (p1)` = fevd_B_p0[, "d4_log_p1"] * 100
)

kable(fevd_table_A, digits = 2,
      caption = "FEVD – Model A: percentage of Δlog(p0) variance explained by own vs agriculture shocks")
kable(fevd_table_B, digits = 2,
      caption = "FEVD – Model B: percentage of Δ4log(p0) variance explained by own vs agriculture shocks")

############################################################
# END OF SCRIPT
############################################################

#New need to rerun gpt explanation

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

library(vars)
roots(var_d1)   # var_model is the output from VAR()
roots(var_d4)   # var_model is the output from VAR()



# ==========================================================
# 1. SETUP DATA & TRANSFORMATIONS
# ==========================================================

# Assuming 'ts_p0' and 'ts_p1' are already loaded in your environment

# --- A. Define Training Cutoff ---
# We use data up to 2023 Q4 for model identification
cut_year <- c(2023, 4)

# --- B. Create LOG Series (Train vs Full) ---
lp0_full <- log(ts_p0)
lp1_full <- log(ts_p1)

# Create TRAINING subsets
lp0_train <- window(lp0_full, end = cut_year)
lp1_train <- window(lp1_full, end = cut_year)

# --- C. Create Differences (on Training Data) ---
# First Differences (1 lag)
dlp0_train <- diff(lp0_train)
dlp1_train <- diff(lp1_train)

# Seasonal Differences (4 lags)
sdlp0_train <- diff(lp0_train, lag = 4)
sdlp1_train <- diff(lp1_train, lag = 4)

# ==========================================================
# 2. DEFINE EXOGENOUS VARIABLES AS TS OBJECTS
# ==========================================================
# Note: We define these for the FULL sample. 
# The testing functions below will automatically cut them 
# to match the Training Data using ts.intersect().

# A. Seasonal Dummies
seas_mat <- seasonaldummy(ts_p0)
seas <- ts(seas_mat, start = start(ts_p0), frequency = frequency(ts_p0))

# B. Pulse Variables (Covid Q2, Covid Q3, Shock 2022 Q2)
tq <- time(ts_p0)
X_pulses_mat <- cbind(
  covid_q2     = as.numeric(tq == 2020 + 1/4),
  covid_q3     = as.numeric(tq == 2020 + 2/4),
  shock_2022q2 = as.numeric(tq == 2022 + 1/4)
)
X_pulses <- ts(X_pulses_mat, start = start(ts_p0), frequency = frequency(ts_p0))

# ==========================================================
# 3. DEFINE ROBUST TESTING FUNCTIONS
# ==========================================================

# --- Robust ADF Function ---
run_adf_robust <- function(y, X_exog, seas_dummies, name="series") {
  cat(paste0("\n=========================================\n"))
  cat(paste0("       ADF TEST for: ", name, "\n"))
  cat(paste0("=========================================\n"))
  
  # 1. ALIGNMENT: 
  # Intersect y (Train), seas (Full), and X (Full).
  # This automatically drops the 'Test' portion of seas/X to match y.
  data_combined <- ts.intersect(y, seas_dummies, X_exog)
  
  # 2. SEPARATE columns back out
  y_aligned <- data_combined[, 1]
  
  # Dynamic column selection
  n_seas <- ncol(seas_dummies)
  seas_aligned <- data_combined[, 2:(1 + n_seas)] 
  X_aligned    <- data_combined[, (2 + n_seas):ncol(data_combined)]
  
  # Define cases
  cases <- list(
    no_tr_no_seas = list(Z = NULL,         deter="none"),
    tr_no_seas    = list(Z = NULL,         deter="trend"),
    seas_no_tr    = list(Z = seas_aligned, deter="none"),
    tr_seas       = list(Z = seas_aligned, deter="trend")
  )
  
  for(nm in names(cases)){
    cat(paste0("\n--- Model: ", nm, " ---\n"))
    Z <- cases[[nm]]$Z
    deter <- cases[[nm]]$deter
    
    # Combine Z (seasonal) and X (pulses)
    if(is.null(Z)) {
      Z_full <- X_aligned
    } else {
      Z_full <- cbind(Z, X_aligned)
    }
    
    # REGRESS OUT EXOGENOUS VARS (Conditional ADF)
    if(!is.null(Z_full)){
      fit    <- lm(y_aligned ~ Z_full)
      y_test <- resid(fit)
    } else {
      y_test <- y_aligned
    }
    
    # Run ADF
    adf <- ur.df(y_test, type = deter, lags = 4, selectlags = "AIC")
    
    # Print clean output
    cat(paste0("Test Statistic: ", round(adf@teststat[1], 4), "\n"))
    cat("Critical Values (tau):\n")
    print(adf@cval[1,]) 
  }
}

# --- Robust KPSS Function ---
run_kpss_robust <- function(y, X_exog, seas_dummies, name="series") {
  cat(paste0("\n========== KPSS for ", name, " ==========\n"))
  
  # 1. ALIGNMENT
  data_combined <- ts.intersect(y, seas_dummies, X_exog)
  y_aligned <- data_combined[, 1]
  
  n_seas <- ncol(seas_dummies)
  seas_aligned <- data_combined[, 2:(1 + n_seas)] 
  X_aligned    <- data_combined[, (2 + n_seas):ncol(data_combined)]
  
  cases <- list(
    no_tr_no_seas = list(Z = NULL,         type="Level"),
    tr_no_seas    = list(Z = NULL,         type="Trend"),
    seas_no_tr    = list(Z = seas_aligned, type="Level"),
    tr_seas       = list(Z = seas_aligned, type="Trend")
  )
  
  for(nm in names(cases)){
    cat(paste0("\n--- ", nm, " ---\n"))
    Z <- cases[[nm]]$Z
    k_type <- cases[[nm]]$type
    
    if(is.null(Z)) {
      Z_full <- X_aligned
    } else {
      Z_full <- cbind(Z, X_aligned)
    }
    
    # Regress out deterministic components
    if(!is.null(Z_full)){
      y_adj <- resid(lm(y_aligned ~ Z_full))
    } else {
      y_adj <- y_aligned
    }
    
    # Run KPSS
    tryCatch({
      k <- kpss.test(y_adj, null = k_type)
      print(k)
    }, error = function(e) {
      cat("KPSS Error (likely convergence):", e$message, "\n")
    })
  }
}

# ==========================================================
# 4. EXECUTE TESTS (ON TRAINING DATA)
# ==========================================================

# ----------------------------------------------------------
# PART A: LOG P0 (GDP) - TRAINING SAMPLE
# ----------------------------------------------------------

# 1. Level (Train)
run_adf_robust(lp0_train, X_pulses, seas, name="log p0 (Level - Train)")
run_kpss_robust(lp0_train, X_pulses, seas, name="log p0 (Level - Train)")

# 2. First Difference (Train)
run_adf_robust(dlp0_train, X_pulses, seas, name="log p0 (First Diff - Train)")
run_kpss_robust(dlp0_train, X_pulses, seas, name="log p0 (First Diff - Train)")

# 3. Seasonal Difference (Train)
run_adf_robust(sdlp0_train, X_pulses, seas, name="log p0 (Seas Diff - Train)")
run_kpss_robust(sdlp0_train, X_pulses, seas, name="log p0 (Seas Diff - Train)")


# ----------------------------------------------------------
# PART B: LOG P1 (Agriculture) - TRAINING SAMPLE
# ----------------------------------------------------------

# 1. Level (Train)
run_adf_robust(lp1_train, X_pulses, seas, name="log p1 (Level - Train)")
run_kpss_robust(lp1_train, X_pulses, seas, name="log p1 (Level - Train)")

# 2. First Difference (Train)
run_adf_robust(dlp1_train, X_pulses, seas, name="log p1 (First Diff - Train)")
run_kpss_robust(dlp1_train, X_pulses, seas, name="log p1 (First Diff - Train)")

# 3. Seasonal Difference (Train)
run_adf_robust(sdlp1_train, X_pulses, seas, name="log p1 (Seas Diff - Train)")
run_kpss_robust(sdlp1_train, X_pulses, seas, name="log p1 (Seas Diff - Train)")


############################################################
# EXTRA. Trend vs Stochastic Trend: Residual Analysis (All Cases)
############################################################

# Load libraries if not already
library(urca)
library(ggplot2)

# --- 1. Prepare Data (Training Set) ---
# Ensure you have these defined from previous steps
# lp0_train, lp1_train (Log Levels)
# seas (Seasonal Dummies - full)

# Create a time index t for the training set
n_train <- length(lp0_train)
t_train <- 1:n_train

# Extract seasonal dummies for train only
# We use 'seas' (ts object) we defined earlier, subsetting it
S_train <- window(seas, end = c(2023, 4)) 
# Convert to matrix for lm()
S_train_mat <- as.matrix(S_train) 
# Drop the last column to avoid collinearity if intercept is present (standard practice)
# or just use it as is if lm() handles singularities. 
# seasonaldummy() usually produces k-1 columns.

# --- 2. Define a Helper Function for Residual Testing ---
test_residuals <- function(y, time_idx, seas_mat, name) {
  cat(paste0("\n>>> Residual Test for: ", name, " <<<\n"))
  
  # Define the 4 Models
  models <- list(
    # 1. Constant Only (No Trend, No Seas)
    none = lm(y ~ 1),
    
    # 2. Linear Trend Only
    trend = lm(y ~ 1 + time_idx),
    
    # 3. Seasonality Only
    seas  = lm(y ~ 1 + seas_mat),
    
    # 4. Trend + Seasonality
    both  = lm(y ~ 1 + time_idx + seas_mat)
  )
  
  results_df <- data.frame(Model=character(), ADF_t_stat=numeric(), Result=character(), stringsAsFactors=FALSE)
  
  for(mod_name in names(models)) {
    # A. Get Residuals
    resids <- resid(models[[mod_name]])
    
    # B. Run ADF on Residuals
    # IMPORTANT: We use type="none" because we have already removed 
    # the deterministic parts (trend/intercept) manually in step A.
    test <- ur.df(resids, type="none", lags=4, selectlags="AIC")
    
    t_stat <- test@teststat[1]
    crit_val <- test@cval[1, 2] # 5% Critical Value
    
    # C. Store Result
    res_text <- ifelse(t_stat < crit_val, "Stationary (Reject H0)", "Non-Stationary (Unit Root)")
    
    cat(sprintf("  %-10s : t=%.4f (5%% CV=%.2f) -> %s\n", 
                mod_name, t_stat, crit_val, res_text))
  }
  
  # Return the residuals of the "Both" model for plotting (usually the most relevant)
  return(resid(models[["both"]]))
}

# --- 3. Execute Tests ---

# A. Run on log p0 (GDP)
res_p0_clean <- test_residuals(lp0_train, t_train, S_train_mat, "log p0 (GDP)")

# B. Run on log p1 (Agriculture)
res_p1_clean <- test_residuals(lp1_train, t_train, S_train_mat, "log p1 (Agriculture)")

# --- 4. Plot the "Cleaned" Residuals (Trend + Seas removed) ---

# Clear any corrupted plots first
if(!is.null(dev.list())) dev.off()

date_train <- time(lp0_train) # Get time index

res_df <- data.frame(
  date = as.numeric(date_train),
  res_p0 = as.numeric(res_p0_clean),
  res_p1 = as.numeric(res_p1_clean)
)

ggplot(res_df, aes(x = date)) +
  geom_line(aes(y = res_p0, colour = "res(p0): Trend+Seas Removed")) +
  geom_line(aes(y = res_p1, colour = "res(p1): Trend+Seas Removed")) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray") +
  labs(title = "Residuals after removing Deterministic Trend & Seasonality",
       subtitle = "Note: If these still show patterns/waves, differencing is preferred.",
       y = "Residual Value", x = "Year", colour = "Series") +
  theme_minimal() +
  theme(legend.position = "bottom")

lp0 <- ts(log_p0, start = c(2015, 1), frequency = 4)

res <- hegy.test(
  lp0,
  deterministic = c(1, 1, 1),        
  lag.method = "AIC",
  maxlag = 4
)

print(res)
summary(res)

lp1 <- ts(log_p1, start = c(2015, 1), frequency = 4)

res <- hegy.test(
  lp1,
  deterministic = c(1, 1, 1),        
  lag.method = "AIC",
  maxlag = 4
)

print(res)
summary(res)

# Test with NO intercept and NO trend
test_none <- ur.df(log_p0, type = "none", lags = 0)
summary(test_none)

# Test with intercept (Drift) but NO trend
test_drift <- ur.df(log_p0, type = "drift", lags = 0)
summary(test_drift)

# Test with intercept and Trend
test_trend <- ur.df(log_p0, type = "trend", lags = 0)
summary(test_trend)

# Test with NO intercept and NO trend
test_none <- ur.df(log_p1, type = "none", lags = 0)
summary(test_none)

# Test with intercept (Drift) but NO trend
test_drift <- ur.df(log_p1, type = "drift", lags = 0)
summary(test_drift)

# Test with intercept and Trend
test_trend <- ur.df(log_p1, type = "trend", lags = 0)
summary(test_trend)


log_p0_stl <- stl(log_p0_tr, s.window = "periodic")
log_p1_stl <- stl(log_p1_tr, s.window = "periodic")

# seasonally adjusted logs
log_p0_sa_tr <- seasadj(log_p0_stl)
log_p1_sa_tr <- seasadj(log_p1_stl)

autoplot(log_p0_stl) + ggtitle("STL decomposition: log_p0_tr")
autoplot(log_p1_stl) + ggtitle("STL decomposition: log_p1_tr")


d1_log_p0_sa_tr <- diff(log_p0_sa_tr)
d1_log_p1_sa_tr <- diff(log_p1_sa_tr)

y_d1_sa <- cbind(d1_log_p0_sa_tr, d1_log_p1_sa_tr)
colnames(y_d1_sa) <- c("d1_log_p0_sa", "d1_log_p1_sa")

# Check basic plots
autoplot(y_d1_sa) +
  ggtitle("First differences of STL-adjusted log GDP and agri GDP") +
  ylab("Δ log (seasonally adjusted)") +
  theme_minimal()

############################################################
# 9C. Model C: VAR in Δ log of seasonally adjusted series
############################################################

# 1-step differences of seasonally adjusted logs
d1_log_p0_sa_tr <- diff(log_p0_sa_tr)
d1_log_p1_sa_tr <- diff(log_p1_sa_tr)

y_d1_sa <- cbind(d1_log_p0_sa_tr, d1_log_p1_sa_tr)
colnames(y_d1_sa) <- c("d1_log_p0_sa", "d1_log_p1_sa")

# Pulses aligned with differenced series (same idea as X_p_tr_d1)
X_pulses_tr      <- X_pulses[train_mask, , drop = FALSE]
X_p_tr_d1        <- X_pulses_tr[-1, , drop = FALSE]   # drop first obs for Δ

# Lag selection without seasonal dummies
lag_sel_d1_sa <- VARselect(y_d1_sa, lag.max = 4,
                           type   = "both",
                           exogen = X_p_tr_d1)
lag_sel_d1_sa$selection
p_d1_sa <- as.integer(lag_sel_d1_sa$selection["AIC(n)"])

var_d1_sa <- VAR(
  y_d1_sa,
  p      = p_d1_sa,
  type   = "both",
  exogen = X_p_tr_d1
)
summary(var_d1_sa)

# Diagnostics
serial.test(var_d1_sa, lags.pt = 8, type = "PT.asymptotic")
arch.test(var_d1_sa, lags.multi = 4)
normality.test(var_d1_sa)

plot(resid(var_d1_sa), main = "VAR(Δ log, STL-adjusted) residuals")



# Run Zivot-Andrews Test
# We use model="both" to allow for a break in intercept AND trend
# lag=NULL allows the function to automatically select lags based on AIC
s1 <- ur.za(log_p0, model = "both", lag = NULL)
s2 <- ur.za(log_p1, model = "both", lag = NULL)
s3 <- ur.za(log_p0, model = "trend", lag = NULL)
s4 <- ur.za(log_p1, model = "trend", lag = NULL)
s5 <- ur.za(log_p0, model = "intercept", lag = NULL)
s6 <- ur.za(log_p1, model = "intercept", lag = NULL)

summary(s1)
summary(s2)
summary(s3)
summary(s4)
summary(s5)
summary(s6)





# Plot the test statistics to see the potential break point
plot(za_test_p0)
plot(za_test_p1)
