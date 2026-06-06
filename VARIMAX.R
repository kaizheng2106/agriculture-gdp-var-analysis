# =========================================================
# --- 0. Load Required Libraries ---
# =========================================================
library(tidyverse)
library(vars)
library(urca)
library(lubridate)
library(tseries) # For adf.test()
library(tsDyn)
library(forecast) # Added for checkresiduals and auto.arima
library(ggplot2) # For advanced plotting

# =========================================================
# --- 1. Load and Wrangle Data ---
# =========================================================
output_directory <- "C:/nottingham 25-26/Math Group Project/CW2 - AGRICULTURE"
file_name <- "original.csv"
full_input_path <- paste0(output_directory, "/", file_name)

# NOTE: Ensure the path is correct or file exists in working directory
if(!file.exists(full_input_path)) {
  cat("Warning: File path not found. Using 'original.csv' in current dir.\n")
  full_input_path <- "original.csv"
}

all_data <- read_csv(full_input_path)

wide_data <- all_data %>%
  filter(series == "abs") %>%
  dplyr::select(-series) %>%
  pivot_wider(names_from = sector, values_from = value) %>%
  mutate(date = ymd(date)) %>%
  arrange(date)

# =========================================================
# --- 2. Create Bivariate Time Series Object (p0, p1) ---
# =========================================================
min_date <- min(wide_data$date)
start_year <- year(min_date)
start_qtr <- quarter(min_date)

ts_data_numeric <- wide_data %>%
  dplyr::select(p0, p1)

data_ts_bivariate <- ts(ts_data_numeric, 
                        start = c(start_year, start_qtr), 
                        frequency = 4)

# =========================================================
# STEP 1: DATA PREPARATION (TRANSFORMATION)
# =========================================================

cat("\n=== STEP 1: DATA PREPARATION ===\n")

# --- MODEL A DATA: Combined Difference (for Raw Currency Changes) ---
# Formula: (1-L)(1-L^4) y_t
p0_comb <- diff(diff(data_ts_bivariate[, "p0"], lag=4), lag=1)
p1_comb <- diff(diff(data_ts_bivariate[, "p1"], lag=4), lag=1)
data_model_A <- na.omit(cbind(p0_comb, p1_comb))
colnames(data_model_A) <- c("p0_CombDiff", "p1_CombDiff")


# --- MODEL B DATA: First Difference of YoY Growth (Acceleration) ---
# YoY Growth: (y_t - y_{t-4}) / y_{t-4}
p0_yoy <- diff(data_ts_bivariate[, "p0"], lag=4) / stats::lag(data_ts_bivariate[, "p0"], k=-4)
p1_yoy <- diff(data_ts_bivariate[, "p1"], lag=4) / stats::lag(data_ts_bivariate[, "p1"], k=-4)
ts_yoy <- na.omit(cbind(p0_yoy, p1_yoy)) # Base YoY Series

# Acceleration Series: First Difference of YoY
p0_accel <- diff(ts_yoy[, "p0_yoy"], lag=1)
p1_accel <- diff(ts_yoy[, "p1_yoy"], lag=1)
data_model_B <- na.omit(cbind(p0_accel, p1_accel))
colnames(data_model_B) <- c("p0_Accel", "p1_Accel")

cat("Model A (Combined Diff) and Model B (Acceleration) data created.\n")
cat("Model B Observations:", nrow(data_model_B), "\n")


# =========================================================
# STEP 5: VAR MODEL A - ESTIMATION
# =========================================================
cat("\n\n================================================================\n")
cat("=== STEP 5: MODEL A - VAR ON COMBINED DIFFERENCED ABSOLUTE DATA ===\n")
cat("================================================================\n")

# 1. Select Lag Order
var_select_A <- VARselect(data_model_A, lag.max = 4, type = "const")
k_A <- min(var_select_A$selection["AIC(n)"], 2)
if(k_A < 1) k_A <- 1
cat("Model A Selected Lag:", k_A, "\n")

# 2. Estimate VAR 
model_A <- VAR(data_model_A, p = k_A, type = "const")
print(summary(model_A))

# 3. Forecast Model A (Just for comparison)
fc_A <- predict(model_A, n.ahead = 4)
plot(fc_A, main="Model A Forecast: Combined Differences")


# =========================================================
# STEP 6: VAR MODEL B - ESTIMATION & EVALUATION
# (The focus model: VAR on Acceleration)
# =========================================================
cat("\n\n================================================================\n")
cat("=== STEP 6: MODEL B - VAR ON FIRST DIFFERENCE OF YoY GROWTH ===\n")
cat("=== VISUAL EVALUATION OF MODEL FIT ===\n")
cat("================================================================\n")

# 1. Select Lag Order
var_select_B <- VARselect(data_model_B, lag.max = 4, type = "const")
k_B <- min(var_select_B$selection["AIC(n)"], 2)
if(k_B < 1) k_B <- 1
cat("Model B Selected Lag:", k_B, "\n")

# --- MODEL FIT EVALUATION SETUP ---
# We use the acceleration series for the test
N_total <- nrow(data_model_B)
N_test <- 8 # Last 8 quarters (2 years) for validation
N_train <- N_total - N_test

ts_train_accel <- data_model_B[1:N_train, ]
test_dates <- time(data_model_B)[(N_train + 1):N_total]

# 2. Retrain VAR on Training Data
model_train <- VAR(ts_train_accel, p = k_B, type = "const")

# 3. Generate 8-Step Forecast (for the known test period)
fc_in_sample <- predict(model_train, n.ahead = N_test)

# 4. Prepare Data for ggplot
df_plot_accel <- data.frame(
  Date = time(data_model_B),
  Actual_p0 = data_model_B[, "p0_Accel"],
  Actual_p1 = data_model_B[, "p1_Accel"],
  Type = factor(c(rep("Train", N_train), rep("Test", N_test)))
)

df_forecast_accel <- data.frame(
  Date = test_dates,
  Forecast_p0 = fc_in_sample$fcst$p0_Accel[, 1],
  Forecast_p1 = fc_in_sample$fcst$p1_Accel[, 1]
)

# --- 5. Plot Comparison (p0 - Total GDP Acceleration) ---
plot_p0_fit <- ggplot(df_plot_accel, aes(x = Date, y = Actual_p0)) +
  # Training Data (Blue)
  geom_line(data = df_plot_accel[df_plot_accel$Type == "Train", ], aes(color = "Training Data"), size = 0.8) +
  # Actual Test Data (Green)
  geom_line(data = df_plot_accel[df_plot_accel$Type == "Test", ], aes(color = "Actual Test Data"), size = 1.2) +
  # Forecast (Dashed Red)
  geom_line(data = df_forecast_accel, aes(x = Date, y = Forecast_p0, color = "Forecast"), 
            linetype = "dashed", size = 1.2) +
  # Mark the split point
  geom_vline(xintercept = time(ts_train_accel)[N_train], linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, color = "gray", linetype = "solid") + # Zero line for acceleration
  scale_x_continuous(breaks = seq(start_year, end(data_model_B)[1], by = 2)) +
  labs(title = "VAR Model B Fit: Total GDP Acceleration (p0)",
       subtitle = paste("Forecast (Dashed Red) compared to Actual data for last", N_test, "quarters."),
       y = "Change in YoY Growth Rate",
       color = "Series") +
  scale_color_manual(values = c("Training Data" = "blue", "Actual Test Data" = "darkgreen", "Forecast" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(plot_p0_fit)


# --- 6. Plot Comparison (p1 - Agriculture Acceleration) ---
plot_p1_fit <- ggplot(df_plot_accel, aes(x = Date, y = Actual_p1)) +
  # Training Data (Blue)
  geom_line(data = df_plot_accel[df_plot_accel$Type == "Train", ], 
            aes(color = "Training Data"), size = 0.8) +
  # Actual Test Data (Green)
  geom_line(data = df_plot_accel[df_plot_accel$Type == "Test", ], 
            aes(color = "Actual Test Data"), size = 1.2) +
  # Forecast (Dashed Red)
  geom_line(data = df_forecast_accel, aes(x = Date, y = Forecast_p1, color = "Forecast"), 
            linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = time(ts_train_accel)[N_train], linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, color = "gray", linetype = "solid") + # Zero line for acceleration
  scale_x_continuous(breaks = seq(start_year, end(data_model_B)[1], by = 2)) +
  labs(title = "VAR Model B Fit: Agriculture GDP Acceleration (p1)",
       subtitle = paste("Forecast (Dashed Red) compared to Actual data for last", N_test, "quarters."),
       y = "Change in YoY Growth Rate",
       color = "Series") +
  scale_color_manual(values = c("Training Data" = "blue", "Actual Test Data" = "darkgreen", "Forecast" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(plot_p1_fit)

cat("\nModel Fit Check Complete: The plots show the VAR model's out-of-sample forecast vs. actual test data.\n")

# =========================================================
# --- 0. Load Required Libraries ---
# =========================================================
library(tidyverse)
library(vars)
library(urca)
library(lubridate)
library(tseries) 
library(tsDyn)
library(forecast) 
library(ggplot2) 

# =========================================================
# --- 1. Load and Wrangle Data ---
# =========================================================
output_directory <- "C:/nottingham 25-26/Math Group Project/CW2 - AGRICULTURE"
file_name <- "original.csv"
full_input_path <- paste0(output_directory, "/", file_name)

# Check file existence
if(!file.exists(full_input_path)) {
  cat("Warning: File path not found. Using 'original.csv' in current dir.\n")
  full_input_path <- "original.csv"
}

all_data <- read_csv(full_input_path)

wide_data <- all_data %>%
  filter(series == "abs") %>%
  dplyr::select(-series) %>%
  pivot_wider(names_from = sector, values_from = value) %>%
  mutate(date = ymd(date)) %>%
  arrange(date)

# =========================================================
# --- 2. Create Bivariate Time Series Object (p0, p1) ---
# =========================================================
min_date <- min(wide_data$date)
start_year <- year(min_date)
start_qtr <- quarter(min_date)

ts_data_numeric <- wide_data %>%
  dplyr::select(p0, p1)

data_ts_bivariate <- ts(ts_data_numeric, 
                        start = c(start_year, start_qtr), 
                        frequency = 4)

# =========================================================
# STEP 1: DATA PREPARATION (YoY Growth & Acceleration)
# =========================================================

cat("\n=== STEP 1: DATA PREPARATION ===\n")

# 1. Calculate YoY Growth (Base Series)
# (y_t - y_{t-4}) / y_{t-4}
p0_yoy <- diff(data_ts_bivariate[, "p0"], lag=4) / stats::lag(data_ts_bivariate[, "p0"], k=-4)
p1_yoy <- diff(data_ts_bivariate[, "p1"], lag=4) / stats::lag(data_ts_bivariate[, "p1"], k=-4)
ts_yoy <- na.omit(cbind(p0_yoy, p1_yoy))

# 2. Calculate Acceleration (First Difference of YoY)
# This is the stationary I(0) data used for the VAR model
p0_accel <- diff(ts_yoy[, "p0_yoy"], lag=1)
p1_accel <- diff(ts_yoy[, "p1_yoy"], lag=1)
data_model_B <- na.omit(cbind(p0_accel, p1_accel))
colnames(data_model_B) <- c("p0_Accel", "p1_Accel")

cat("Model B (Acceleration) data created. Observations:", nrow(data_model_B), "\n")


# =========================================================
# STEP 6: VAR MODEL B - ESTIMATION
# =========================================================
cat("\n=== STEP 6: ESTIMATING VAR MODEL B ===\n")

# 1. Select Lag Order 
var_select_B <- VARselect(data_model_B, lag.max = 4, type = "const")
k_B <- min(var_select_B$selection["AIC(n)"], 2)
if(k_B < 1) k_B <- 1
cat("Selected Lag Order (AIC):", k_B, "\n")

# 2. Define Train/Test Split (Last 6 Quarters: 2024 Q1 - 2025 Q2)
N_total <- nrow(data_model_B)
N_test <- 6 
N_train <- N_total - N_test

ts_train_accel <- data_model_B[1:N_train, ]
ts_test_accel <- data_model_B[(N_train + 1):N_total, ]
test_dates <- time(data_model_B)[(N_train + 1):N_total]

# 3. Train VAR on Training Set
model_train <- VAR(ts_train_accel, p = k_B, type = "const")

# 4. Generate Forecast
fc_in_sample <- predict(model_train, n.ahead = N_test)

# 5. Calculate Accuracy (RMSE)
actual_p0 <- ts_test_accel[, "p0_Accel"]
forecast_p0 <- fc_in_sample$fcst$p0_Accel[, 1]
rmse_p0 <- sqrt(mean((actual_p0 - forecast_p0)^2))

actual_p1 <- ts_test_accel[, "p1_Accel"]
forecast_p1 <- fc_in_sample$fcst$p1_Accel[, 1]
rmse_p1 <- sqrt(mean((actual_p1 - forecast_p1)^2))

cat(paste("\nAccuracy Metric (RMSE):\n"))
cat(paste("  Total GDP (p0):", round(rmse_p0, 4), "\n"))
cat(paste("  Agriculture (p1):", round(rmse_p1, 4), "\n"))


# =========================================================
# STEP 7: VISUALIZATION (Pattern Match)
# =========================================================
cat("\n=== STEP 7: GENERATING PLOTS ===\n")

# Define the custom plotting function matching your requested pattern
plot_fc_var <- function(series_name, actual_series, forecast_values, rmse_value) {
  
  # Prepare data subsets
  train_series <- actual_series[1:N_train]
  test_series <- actual_series[(N_train + 1):N_total]
  
  # Construct main plotting dataframe
  df_plot <- data.frame(
    Date = c(time(ts_train_accel), test_dates, test_dates),
    Value = c(as.numeric(train_series),
              as.numeric(forecast_values),
              as.numeric(test_series)),
    Type = c(rep("Train (2016–2023)", length(train_series)),
             rep("Forecast (2024–2025)", length(forecast_values)),
             rep("Actual Test (2024–2025)", length(test_series)))
  )
  
  # Create a connecting point for the forecast line so it starts from the last training point
  last_train_point <- data.frame(
    Date = time(ts_train_accel)[N_train],
    Value = as.numeric(ts_train_accel)[N_train, which(colnames(data_model_B) == colnames(actual_series))], # Ensure correct column selection isn't needed here as passed vector is 1D, but safe logic
    Type = "Forecast (2024–2025)"
  )
  # Fix: actual_series is a vector here, so just take the last element
  last_train_point$Value <- as.numeric(train_series)[length(train_series)]
  
  df_forecast_continuous <- bind_rows(last_train_point, df_plot[df_plot$Type == "Forecast (2024–2025)", ])
  
  # Generate Plot
  p <- ggplot(df_plot, aes(x=Date, y=Value, color=Type)) +
    # 1. Train Data (Grey)
    geom_line(data=df_plot[df_plot$Type == "Train (2016–2023)", ], 
              aes(x=Date, y=Value), size=1.1, color="gray40") +
    
    # 2. Actual Test Data (Black)
    geom_line(data=df_plot[df_plot$Type == "Actual Test (2024–2025)", ], 
              aes(x=Date, y=Value), size=1.1, color="black") +
    
    # 3. Forecast Data (Blue)
    geom_line(data=df_forecast_continuous, 
              aes(x=Date, y=Value), size=1.1, color="blue") +
    
    # Define Legend Colors
    scale_color_manual(values=c(
      "Train (2016–2023)"="gray40",
      "Forecast (2024–2025)"="blue",
      "Actual Test (2024–2025)"="black"
    )) +
    
    # Add Reference Lines
    geom_vline(xintercept = time(ts_train_accel)[N_train], linetype = "dotted", color = "red") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dotted") +
    
    # Labels and Theme
    labs(
      title=paste0("Forecast vs Actual: ", series_name),
      subtitle=paste("Model: VAR(", k_B, ") on Acceleration. RMSE:", round(rmse_value, 4)),
      y="Change in YoY Growth Rate",
      x="Year",
      color="Series"
    ) +
    theme_minimal(base_size=14) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold"),
      legend.position="bottom"
    )
  
  return(p)
}

# --- Generate and Print Plots ---

# 1. Plot for Total GDP (p0)
plot_p0 <- plot_fc_var(
  series_name = "Total GDP (p0)",
  actual_series = data_model_B[, "p0_Accel"],
  forecast_values = forecast_p0,
  rmse_value = rmse_p0
)
print(plot_p0)

# 2. Plot for Agriculture GDP (p1)
plot_p1 <- plot_fc_var(
  series_name = "Agriculture GDP (p1)",
  actual_series = data_model_B[, "p1_Accel"],
  forecast_values = forecast_p1,
  rmse_value = rmse_p1
)
print(plot_p1)

cat("\nSuccess: Forecast vs. Actual plots generated.\n")