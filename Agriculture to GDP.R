# ================================================
# Malaysia GDP: Contribution of Agriculture Sector 
# ================================================

# ---- 1. Load packages ----
library(tidyverse)
library(tseries)
library(forecast)
library(TSA)
library(ggplot2)
library(gridExtra)
library(tsoutliers)

# ---- 2. Read and prepare data ----
output_directory <- "C:/nottingham 25-26/Math Group Project/CW2 - AGRICULTURE"
file_name <- "original.csv"
full_input_path <- paste0(output_directory, "/", file_name)

gdp <- read_csv(full_input_path) 
gdp$date <- as.Date(gdp$date)

# --- A. Ratio Data Preparation ---
df_abs <- gdp %>%
  filter(series == "abs", sector %in% c("p0", "p1")) %>%
  pivot_wider(names_from = sector, values_from = value) %>%
  arrange(date) %>%
  mutate(ratio = p1 / p0)

y_ts <- ts(df_abs$ratio, start = c(2015, 1), frequency = 4)
n_obs_abs <- length(y_ts)


# =========================================================
# STEP 3: STATIONARITY AND VISUAL DIAGNOSTICS
# =========================================================

cat("\n\n=== STEP 3: STATIONARITY JUSTIFICATION AND VISUAL DIAGNOSTICS ===\n")

# --- 3.1 Plot RAW Data (Level and Share) ---
p_total <- ggplot(df_abs, aes(x = date, y = p0)) + geom_line(color = "blue", linewidth = 1) + labs(title = "Total GDP (p0) Level", y = "Nominal GDP", x = "Quarter") + theme_minimal(base_size = 14)
p_agri <- ggplot(df_abs, aes(x = date, y = p1)) + geom_line(color = "darkgreen", linewidth = 1) + labs(title = "Agriculture GDP (p1) Level", y = "Nominal GDP", x = "Quarter") + theme_minimal(base_size = 14)
grid.arrange(p_total, p_agri, ncol = 1)

p_ratio <- ggplot(df_abs, aes(x = date, y = ratio)) + geom_line(color = "purple", linewidth = 1) + labs(title = "Agriculture Share of Total GDP (p1/p0)", y = "Ratio", x = "Quarter") + theme_minimal(base_size = 14)
print(p_ratio)


# --- 3.2 STEP-BY-STEP ADF Tests (Full Justification) ---
cat("\n--- 3.2 FULL ADF TEST OUTPUTS (Ratio Only) ---\n")

# 1. RAW Ratio (P1/P0) - Must be Non-stationary
cat("\n[1. RAW RATIO]\n"); print(adf.test(df_abs$ratio)) 

# 2. Ratio Differencing Steps (Justifying d=1, D=1)
ratio_diff_s <- diff(df_abs$ratio, lag = 4)
cat("\n[2. SEASONAL DIFFERENCE (D=1)]\n"); print(adf.test(na.omit(ratio_diff_s)))

ratio_combined_diff <- diff(ratio_diff_s, lag = 1)
cat("\n[3. COMBINED DIFFERENCE (d=1, D=1)]\n"); print(adf.test(na.omit(ratio_combined_diff)) )
f1 <- na.omit(ratio_combined_diff) # Final stationary ratio series

# --- 3.3 ACF, PACF, EACF (on the stationary ratio series) ----
cat("\n--- 3.3 ACF, PACF, and EACF on Stationary Ratio (d=1, D=1) ---\n")
par(mfrow = c(1, 2))
acf(f1, main = "ACF of Stationary Ratio (d=1, D=1)", lag.max = 40)
pacf(f1, main = "PACF of Stationary Ratio (d=1, D=1)", lag.max = 40)
par(mfrow = c(1, 1)) # Reset plot window
eacf(f1, ar.max = 8, ma.max = 9)


# --- 3.4 Plotting Stationary Series (CONFIRMATION) ---
par(mfrow = c(1, 1))
ts.plot(f1, main = "Combined Difference of P1/P0 Ratio (Stationary)", ylab = "Delta(Delta_4 Ratio)", col = "purple")
par(mfrow = c(1, 1)) 


# =========================================================
# STEP 4: XREG CREATION (Outlier Intervention Variables)
# =========================================================

cat("\n\n=== STEP 4: OUTLIER DETECTION AND ANALYSIS ===\n")

# --- 4.1 Detect Outliers on Original Ratio Series ---
outliers_share <- tsoutliers(y_ts)

cat("\n--- Outlier Detection Results ---\n")
print(outliers_share)

# Display outlier details
if (length(outliers_share$index) > 0) {
  outlier_df <- data.frame(
    Index = outliers_share$index,
    Date = time(y_ts)[outliers_share$index],
    Type = outliers_share$type,
    Coefficient = outliers_share$coefhat
  )
  cat("\nDetected Outliers:\n")
  print(outlier_df)
  
  # Plot outliers on the time series
  plot(y_ts, main = "Agriculture Share with Detected Outliers", 
       ylab = "Ratio", xlab = "Time", col = "blue", lwd = 2)
  points(time(y_ts)[outliers_share$index], 
         y_ts[outliers_share$index], 
         col = "red", pch = 19, cex = 2)
  legend("topright", legend = c("Original Series", "Outliers"), 
         col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 19))
  
} else {
  cat("\nNo outliers detected.\n")
}

# --- 4.2 Additional Outlier Detection Methods (for comparison) ---

# Z-score method on stationary series (f1)
z_scores <- abs(scale(f1))
z_outliers <- which(z_scores > 3)

cat("\n--- Z-Score Method (threshold = 3) on Stationary Series ---\n")
if (length(z_outliers) > 0) {
  cat("Outliers detected at indices:", z_outliers, "\n")
  cat("Corresponding dates:", time(y_ts)[z_outliers + 5], "\n")  # Adjust for differencing
} else {
  cat("No outliers detected using Z-score method.\n")
}

# IQR method on stationary series
Q1 <- quantile(f1, 0.25)
Q3 <- quantile(f1, 0.75)
IQR_val <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val
iqr_outliers <- which(f1 < lower_bound | f1 > upper_bound)

cat("\n--- IQR Method on Stationary Series ---\n")
cat("Lower bound:", lower_bound, "Upper bound:", upper_bound, "\n")
if (length(iqr_outliers) > 0) {
  cat("Outliers detected at indices:", iqr_outliers, "\n")
} else {
  cat("No outliers detected using IQR method.\n")
}

# --- 4.3 Visual Diagnostics ---
par(mfrow = c(2, 2))

# Boxplot of stationary series
boxplot(f1, main = "Boxplot of Stationary Series", 
        ylab = "Combined Differenced Ratio", col = "lightblue")

# Histogram with normal curve overlay
hist(f1, probability = TRUE, main = "Histogram of Stationary Series",
     xlab = "Value", col = "lightgreen", breaks = 15)
curve(dnorm(x, mean = mean(f1), sd = sd(f1)), add = TRUE, col = "red", lwd = 2)

# Q-Q plot
qqnorm(f1, main = "Q-Q Plot of Stationary Series")
qqline(f1, col = "red", lwd = 2)

# Time series with control limits
plot(f1, main = "Stationary Series with Control Limits",
     ylab = "Value", type = "l")
abline(h = mean(f1) + 3*sd(f1), col = "red", lty = 2, lwd = 2)
abline(h = mean(f1) - 3*sd(f1), col = "red", lty = 2, lwd = 2)
abline(h = mean(f1), col = "blue", lty = 1)

par(mfrow = c(1, 1))

# --- 4.4 Create Outlier Regressors (Your existing code continues) ---
outlier_indices <- outliers_share$index
n_outliers <- length(outlier_indices)
outlier_regressors_share <- NULL 

if (n_outliers > 0) {
  outlier_regressors_share <- matrix(0, nrow = n_obs_abs, ncol = n_outliers)
  for (i in 1:n_outliers) {
    k <- outlier_indices[i] 
    if (k <= n_obs_abs) {
      outlier_regressors_share[k, i] <- 1
    }
  }
  colnames(outlier_regressors_share) <- paste0("outlier_", outlier_indices)
  
  cat("\n--- Outlier Regressors Created ---\n")
  cat("Number of outlier regressors:", ncol(outlier_regressors_share), "\n")
} else {
  cat("\nNo outlier regressors needed.\n")
}

# =========================================================
# ADDITIONAL ANALYSIS: Outlier Detection on Raw p0 and p1
# =========================================================

cat("\n\n=== OUTLIER DETECTION ON RAW GDP SERIES ===\n")

# Create time series for p0 and p1
p0_ts <- ts(df_abs$p0, start = c(2015, 1), frequency = 4)
p1_ts <- ts(df_abs$p1, start = c(2015, 1), frequency = 4)

# Detect outliers in p0 (Total GDP)
cat("\n--- Outliers in Total GDP (p0) ---\n")
outliers_p0 <- tsoutliers(p0_ts)
print(outliers_p0)
if (length(outliers_p0$index) > 0) {
  cat("Number of outliers in p0:", length(outliers_p0$index), "\n")
  cat("Outlier indices:", outliers_p0$index, "\n")
  cat("Outlier types:", outliers_p0$type, "\n")
}

# Detect outliers in p1 (Agriculture GDP)
cat("\n--- Outliers in Agriculture GDP (p1) ---\n")
outliers_p1 <- tsoutliers(p1_ts)
print(outliers_p1)
if (length(outliers_p1$index) > 0) {
  cat("Number of outliers in p1:", length(outliers_p1$index), "\n")
  cat("Outlier indices:", outliers_p1$index, "\n")
  cat("Outlier types:", outliers_p1$type, "\n")
}

# Visual comparison
par(mfrow = c(3, 1))
plot(p0_ts, main = "Total GDP (p0)", ylab = "Nominal GDP", col = "blue", lwd = 2)
if (length(outliers_p0$index) > 0) {
  points(time(p0_ts)[outliers_p0$index], p0_ts[outliers_p0$index], 
         col = "red", pch = 19, cex = 2)
}

plot(p1_ts, main = "Agriculture GDP (p1)", ylab = "Nominal GDP", col = "darkgreen", lwd = 2)
if (length(outliers_p1$index) > 0) {
  points(time(p1_ts)[outliers_p1$index], p1_ts[outliers_p1$index], 
         col = "red", pch = 19, cex = 2)
}

plot(y_ts, main = "Agriculture Share (p1/p0)", ylab = "Ratio", col = "purple", lwd = 2)
par(mfrow = c(1, 1))

# Summary comparison
cat("\n--- SUMMARY ---\n")
cat("Outliers in p0 (Total GDP):", length(outliers_p0$index), "\n")
cat("Outliers in p1 (Agriculture GDP):", length(outliers_p1$index), "\n")
cat("Outliers in p1/p0 (Ratio):", length(outliers_share$index), "\n")


# =========================================================
# STEP 5: MODEL FITTING (Candidates on Full Data)
# =========================================================

cat("\n\n=== STEP 5: MODEL FITTING (Share Models on Full Data) ===\n")

# --- Candidate 1: Raw Ratio (Auto d, D) ---
fit_share_raw <- auto.arima(y = y_ts, xreg = outlier_regressors_share, seasonal = TRUE, stepwise = FALSE, trace = FALSE)
cat("\n--- Model 1A (Share/Auto): Best Model Found ---\n"); print(summary(fit_share_raw))


# --- Candidate 2: Raw Ratio (Fixed d=1, D=1) ---
fit_share_fixed <- auto.arima(y = y_ts, d=1, D=1, xreg = outlier_regressors_share, seasonal = TRUE, stepwise = FALSE, trace = FALSE)
cat("\n--- Model 1B (Share/Fixed d=1, D=1): Model Forced to Manual Differencing ---\n"); print(summary(fit_share_fixed))


# =========================================================
# STEP 6: COMPARISON, DIAGNOSTICS, AND FORECASTING
# =========================================================

# --- 6.1 Helper Function (Fixed) ---
diagnose_model <- function(model, model_name, test_data = NULL, xreg_future = NULL) {
  mse <- NA
  aic_val <- NA
  lb_p <- NA
  
  if (!is.null(test_data) && length(test_data) > 0) {
    fc <- forecast(model, h = length(test_data), xreg = xreg_future)
    mse <- mean((test_data - fc$mean)^2) 
  }
  
  is_arima_model <- inherits(model, "Arima")
  
  if (is_arima_model) {
    aic_val <- AIC(model)
    lb <- Box.test(residuals(model), lag = 20, type = "Ljung-Box", fitdf = length(model$coef))
    lb_p <- lb$p.value
  } else {
    lb <- Box.test(residuals(model), lag = 20, type = "Ljung-Box")
    lb_p <- lb$p.value
  }
  
  return(data.frame(
    Model = model_name,
    AIC = aic_val,
    LjungBox_p = lb_p,
    MSE_Test = mse
  ))
}

# --- 6.2 Train–Test Split and Refitting (for Validation) ---
train_end <- c(2023, 4)
test_start <- c(2024, 1)

train_ratio <- window(y_ts, end = train_end)
test_ratio <- window(y_ts, start = test_start)
train_len_abs <- length(train_ratio)

# Prepare XREG splits
if (!is.null(outlier_regressors_share)) {
  train_X_reg_share <- outlier_regressors_share[1:train_len_abs, ]
  test_X_reg_share <- matrix(0, nrow = length(test_ratio), ncol = ncol(outlier_regressors_share))
  colnames(test_X_reg_share) <- colnames(outlier_regressors_share)
} else {
  train_X_reg_share <- NULL
  test_X_reg_share <- NULL
}

# Refit on Training Data (Essential step for comparison)
fit_train_share_raw <- auto.arima(train_ratio, xreg = train_X_reg_share, seasonal = TRUE, stepwise = FALSE, trace=FALSE)
fit_train_share_fixed <- auto.arima(train_ratio, d=1, D=1, xreg = train_X_reg_share, seasonal = TRUE, stepwise = FALSE, trace=FALSE)
fit_train_snaive_share <- snaive(train_ratio) # Benchmark


# --- 6.3 Final Comparison Table ---
cat("\n\n=== FINAL COMPARISON OF CANDIDATE MODELS ===\n")
cat("NOTE: LjungBox_p > 0.05 is REQUIRED for model validity. Lowest MSE is best accuracy.\n")

comp_table_final <- rbind(
  diagnose_model(fit_train_share_raw, "1A: SARIMAX (Share/Auto d,D)", test_ratio, test_X_reg_share),
  diagnose_model(fit_train_share_fixed, "1B: SARIMAX (Share/Fixed d=1,D=1)", test_ratio, test_X_reg_share),
  diagnose_model(fit_train_snaive_share, "1C: SNAIVE (Share Benchmark)", test_ratio)
)
print(comp_table_final)


# --- 6.4 Full Residual Diagnostics (For Validation and Final Report) ---
cat("\n\n=== RESIDUAL DIAGNOSTICS FOR BEST MODELS ===\n")

best_share_model <- fit_train_share_fixed
cat("\n--- Model 1B: SARIMAX (Structural Share/Fixed d=1,D=1) Residuals ---\n")
res_share <- residuals(best_share_model)
par(mfrow = c(2, 2))
ts.plot(res_share, main = "Share Residuals over Time", ylab="Residuals")
acf(res_share, main = "Share Residual ACF", lag.max = 20)
pacf(res_share, main = "Share Residual PACF", lag.max = 20)
hist(res_share, main = "Share Residual Histogram", col = "lightblue")
par(mfrow = c(1, 1)) 
checkresiduals(best_share_model)


# --- 6.5 Forecast Visualization (8 Quarters) ---
h <- 8  

cat("\n\n=== FINAL FORECAST VISUALIZATIONS ===\n")

future_Xreg_share <- matrix(0, nrow = h, ncol = ncol(test_X_reg_share))
colnames(future_Xreg_share) <- colnames(test_X_reg_share)

# Forecast for Structural Share (Using the statistically justified fixed model)
safe_forecast <- function(model, h, xreg_future) {
  if (is.null(xreg_future)) {
    return(forecast(model, h = h))
  } else {
    return(forecast(model, h = h, xreg = xreg_future))
  }
}

f_share <- safe_forecast(best_share_model, h, future_Xreg_share)
plot_share <- autolayer(f_share) + autolayer(y_ts, series = "Observed", color = "black") + 
  labs(title = "Forecast: Agriculture Share of GDP (p1/p0)", y = "Ratio", x = "Quarter") + theme_minimal(base_size = 14)
print(plot_share)

fit_share_fixed <- auto.arima(y = y_ts, d=1, D=1, seasonal = TRUE, stepwise = FALSE, trace = FALSE)
print(fit_share_fixed)








# ============================
# FIT SARIMA(0,1,0)(0,1,0)[4] WITH XREG
# ============================

fit_sarima010_010_xreg <- Arima(
  y_ts,
  order = c(0,1,0),
  seasonal = list(order = c(0,1,0), period = 4),
  xreg = outlier_regressors_share,
  include.drift = FALSE,   # SARIMA(0,1,0)(0,1,0) has NO drift
  method = "ML"
)

summary(fit_sarima010_010_xreg)
checkresiduals(fit_sarima010_010_xreg)

# ============================
# CREATE FUTURE XREG (all zero)
# ============================

h <- 8  # forecast horizon

if (!is.null(outlier_regressors_share)) {
  future_xreg <- matrix(0, nrow = h, ncol = ncol(outlier_regressors_share))
  colnames(future_xreg) <- colnames(outlier_regressors_share)
} else {
  future_xreg <- NULL
}

# ============================
# FORECASTING
# ============================

fc_sarima010_010_xreg <- forecast(
  fit_sarima010_010_xreg,
  h = h,
  xreg = future_xreg
)

print(fc_sarima010_010_xreg)

autoplot(fc_sarima010_010_xreg) +
  autolayer(y_ts, series = "Observed") +
  labs(
    title = "Forecast: SARIMA(0,1,0)(0,1,0)[4] with Outlier XREG",
    y = "Agriculture Share (p1/p0)",
    x = "Quarter"
  ) +
  theme_minimal(base_size = 14)

# ============================
# TRAIN / TEST SPLIT (2015–2023 train, 2024–2025 test)
# ============================

train_end <- c(2023, 4)
test_start <- c(2024, 1)

train_ratio <- window(y_ts, end = train_end)
test_ratio  <- window(y_ts, start = test_start)

train_len <- length(train_ratio)
test_len  <- length(test_ratio)

# Split XREG
if (!is.null(outlier_regressors_share)) {
  train_xreg <- outlier_regressors_share[1:train_len, ]
  test_xreg  <- outlier_regressors_share[(train_len+1):(train_len + test_len), ]
} else {
  train_xreg <- NULL
  test_xreg  <- NULL
}

# ============================
# PLOT: TRAIN vs FORECAST vs ACTUAL (CLEAR LEGEND)
# ============================

# Build forecast object with explicit ts labels
fc_test_plot <- fc_test
fc_test_plot$mean <- ts(fc_test$mean, start = test_start, frequency = 4)

# Create dataframe for manual ggplot
df_plot <- data.frame(
  Date = c(time(train_ratio), time(fc_test_plot$mean), time(test_ratio)),
  Value = c(as.numeric(train_ratio), as.numeric(fc_test_plot$mean), as.numeric(test_ratio)),
  Type = c(
    rep("Train (2015–2023)", length(train_ratio)),
    rep("Forecast (2024–2025)", length(fc_test_plot$mean)),
    rep("Actual Test (2024–2025)", length(test_ratio))
  )
)

# Plot
ggplot(df_plot, aes(x = Date, y = Value, color = Type)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "Train (2015–2023)" = "gray40",
      "Forecast (2024–2025)" = "blue",
      "Actual Test (2024–2025)" = "black"
    )
  ) +
  labs(
    title = "SARIMA(0,1,0)(0,1,0)[4] with XREG\nForecast vs Actual",
    x = "Year",
    y = "Agriculture Share (p1/p0)",
    color = "Series"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# ------------------------------
# 0. Prepare Train/Test Split
# ------------------------------
train_end <- c(2023, 4)
test_start <- c(2024, 1)

train_ratio <- window(y_ts, end = train_end)
test_ratio  <- window(y_ts, start = test_start)

# ------------------------------
# 1. Fit FOUR SARIMA MODELS (NO XREG)
# ------------------------------

m1 <- Arima(train_ratio, order=c(0,1,0),
            seasonal=list(order=c(1,1,1), period=4))

m2 <- Arima(train_ratio, order=c(0,1,1),
            seasonal=list(order=c(1,1,1), period=4))

m3 <- Arima(train_ratio, order=c(0,1,2),
            seasonal=list(order=c(1,1,1), period=4))

m4 <- Arima(train_ratio, order=c(1,1,1),
            seasonal=list(order=c(1,1,1), period=4))

models <- list(m1=m1, m2=m2, m3=m3, m4=m4)

# ------------------------------
# 2. Forecast Test Set (2024–2025)
# ------------------------------

fc <- lapply(models, forecast, h=length(test_ratio))

# ------------------------------
# 3. Function For Train-Forecast-Actual PLOT
# ------------------------------

plot_fc <- function(model_name, model_fc) {
  
  df_plot <- data.frame(
    Date = c(time(train_ratio), time(model_fc$mean), time(test_ratio)),
    Value = c(as.numeric(train_ratio),
              as.numeric(model_fc$mean),
              as.numeric(test_ratio)),
    Type = c(rep("Train (2015–2023)", length(train_ratio)),
             rep("Forecast (2024–2025)", length(model_fc$mean)),
             rep("Actual Test (2024–2025)", length(test_ratio)))
  )
  
  ggplot(df_plot, aes(x=Date, y=Value, color=Type)) +
    geom_line(size=1.1) +
    scale_color_manual(values=c(
      "Train (2015–2023)"="gray40",
      "Forecast (2024–2025)"="blue",
      "Actual Test (2024–2025)"="black"
    )) +
    labs(
      title=paste0("Forecast vs Actual: ", model_name),
      y="Agriculture Share (p1/p0)",
      x="Year",
      color="Series"
    ) +
    theme_minimal(base_size=14) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold"),
      legend.position="bottom"
    )
}

# ------------------------------
# 4. PRINT PLOTS FOR ALL FOUR MODELS
# ------------------------------

plot_fc("SARIMA(0,1,0)(1,1,1)[4]", fc$m1)
plot_fc("SARIMA(0,1,1)(1,1,1)[4]", fc$m2)
plot_fc("SARIMA(0,1,2)(1,1,1)[4]", fc$m3)
plot_fc("SARIMA(1,1,1)(1,1,1)[4]", fc$m4)

# ------------------------------
# 5. MODEL COMPARISON TABLE (AIC/BIC + Test MSE)
# ------------------------------

compare <- data.frame(
  Model = c("SARIMA(0,1,0)(1,1,1)[4]",
            "SARIMA(0,1,1)(1,1,1)[4]",
            "SARIMA(0,1,2)(1,1,1)[4]",
            "SARIMA(1,1,1)(1,1,1)[4]"),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  Test_MSE = sapply(1:4, function(i){
    mean((test_ratio - fc[[i]]$mean)^2)
  })
)

print(compare)

# --- Step 1: Fit SARIMAX on TRAIN data ---
fit_train_sarima010_010_xreg <- Arima(
  train_ratio,
  order = c(0,1,0),
  seasonal = list(order = c(0,1,0), period = 4),
  xreg = train_xreg,
  include.drift = FALSE
)

# --- Step 2: Forecast for TEST period (2024–2025) ---
fc_010010_xreg <- forecast(
  fit_train_sarima010_010_xreg,
  h = length(test_ratio),
  xreg = test_xreg      # IMPORTANT!!
)

# --- Step 3: Accuracy ---

accuracy(fc$m1, test_ratio)
# accuracy(fc$m2, test_ratio)
# accuracy(fc$m3, test_ratio)
# accuracy(fc$m4, test_ratio)
accuracy(fc_010010_xreg, test_ratio)


# ============================
# SETUP: TRAIN/TEST SPLIT
# ============================

train_end <- c(2023, 4)
test_start <- c(2024, 1)

train_ratio <- window(y_ts, end = train_end)
test_ratio  <- window(y_ts, start = test_start)

train_len <- length(train_ratio)
test_len  <- length(test_ratio)

# Split XREG if available
if (!is.null(outlier_regressors_share)) {
  train_xreg <- outlier_regressors_share[1:train_len, ]
  test_xreg  <- outlier_regressors_share[(train_len+1):(train_len + test_len), ]
} else {
  train_xreg <- NULL
  test_xreg  <- NULL
}

# ============================
# FIT SARIMA(0,1,0)(0,1,0)[4] WITH XREG ON TRAIN DATA
# ============================

fit_sarima <- Arima(
  train_ratio,
  order = c(0, 1, 0),
  seasonal = list(order = c(0, 1, 0), period = 4),
  xreg = train_xreg,
  include.drift = FALSE,
  method = "ML"
)

summary(fit_sarima)
checkresiduals(fit_sarima)

# ============================
# FORECAST 2024–2025 (TEST PERIOD)
# ============================

fc_2024_2025 <- forecast(
  fit_sarima,
  h = test_len,
  xreg = test_xreg
)

print(fc_2024_2025)

# ============================
# PLOT: TRAIN vs FORECAST vs ACTUAL
# ============================

df_plot <- data.frame(
  Date = c(time(train_ratio), time(fc_2024_2025$mean), time(test_ratio)),
  Value = c(
    as.numeric(train_ratio),
    as.numeric(fc_2024_2025$mean),
    as.numeric(test_ratio)
  ),
  Type = c(
    rep("Train (2015–2023)", length(train_ratio)),
    rep("Forecast (2024–2025)", length(fc_2024_2025$mean)),
    rep("Actual (2024–2025)", length(test_ratio))
  )
)

ggplot(df_plot, aes(x = Date, y = Value, color = Type)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "Train (2015–2023)" = "gray40",
      "Forecast (2024–2025)" = "blue",
      "Actual (2024–2025)" = "black"
    )
  ) +
  labs(
    title = "SARIMA(0,1,0)(0,1,0)[4] with XREG\nForecast vs Actual (2024–2025)",
    x = "Year",
    y = "Agriculture Share (p1/p0)",
    color = "Series"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# ============================
# ACCURACY METRICS
# ============================

accuracy(fc_2024_2025, test_ratio)