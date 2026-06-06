# Install necessary packages if you haven't already
# install.packages(c("readr", "dplyr", "vars", "urca"))

# Load libraries
library(readr)
library(dplyr)
library(vars)
library(urca)

# --- 0. Data Preparation ---

# 1. Load the data
df <- read_csv("original.csv")

# 2. Filter for 'abs' series, pivot, and convert to time series (ts) object
# Assuming quarterly data, starting Q1 2015
data_abs_wide <- df %>%
  filter(series == "abs", sector %in% c("p0", "p1")) %>%
  select(date, sector, value) %>%
  tidyr::pivot_wider(names_from = sector, values_from = value) %>%
  mutate(date = as.Date(date)) %>%
  select(date, p0, p1) %>%
  na.omit()

# Convert to ts object for VAR/VECM packages (start year, frequency=4 for quarterly)
# Adjust start parameter if your first date is not 2015-01-01
start_year <- as.numeric(format(min(data_abs_wide$date), "%Y"))
start_quarter <- as.numeric(format(min(data_abs_wide$date), "%m")) / 3
ts_data_levels <- ts(data_abs_wide[, c("p0", "p1")], 
                     start = c(start_year, start_quarter), 
                     frequency = 4)

# 3. Create Exogenous Dummies (Seasonal and Pulse)
# Seasonal Dummies (Q1 is omitted, so we need Q2, Q3, Q4)
Q2 <- as.numeric(cycle(ts_data_levels) == 2)
Q3 <- as.numeric(cycle(ts_data_levels) == 3)
Q4 <- as.numeric(cycle(ts_data_levels) == 4)

# Pulse Dummies (assuming Q2 2020 and Q3 2022)
# Date index: ts_data_levels[i,] corresponds to the i-th row of the original data.
# We need to find the index for the shock dates.
shock_dates <- c("2020-04-01", "2022-07-01")
shock_indices <- which(data_abs_wide$date %in% as.Date(shock_dates))
D_COVID <- rep(0, nrow(ts_data_levels))
D_2022_SHOCK <- rep(0, nrow(ts_data_levels))

if (length(shock_indices) > 0) {
  D_COVID[shock_indices[1]] <- 1 # Assuming first shock date is COVID
  D_2022_SHOCK[shock_indices[2]] <- 1 # Assuming second shock date is 2022
}

exog_vars <- data.frame(Q2, Q3, Q4, D_COVID, D_2022_SHOCK)


# ----------------------------------------------------------------------------------
# --- Step 1: Check for Stationarity (Combined Difference) ---
# We already confirmed I(1,1) in the previous run, but this confirms I(0) for the
# stationary modeling base.

cat("\n--- Step 1: ADF Test on Combined Difference (I(0) Check) ---\n")

# Calculate the combined difference: Delta_1 * Delta_4 * y_t
ts_data_combined_diff <- diff(ts_data_levels, lag = 1) %>% 
  diff(lag = 4) %>%
  na.omit()

# ADF test on combined difference (regression="none" for I(0) series)
# We test each variable individually.
for (var in c("p0", "p1")) {
  adf_test <- ur.df(ts_data_combined_diff[, var], type = "none", lags = 1)
  cat(paste0(var, " (Combined Diff): \n"))
  print(summary(adf_test))
  # Interpretation: If Test Statistic < Critical Value, Reject H0 (Unit Root)
}

# ----------------------------------------------------------------------------------
# --- Step 2: Select the Optimal Lag Length (p) ---
# Use VAR on LEVELS data with DUMMIES up to max lag 4 (Lecturer's Advice)

# Must align the exogenous variables with the time series length (some rows are removed by diff in later steps)
exog_vars_aligned <- exog_vars[1:nrow(ts_data_levels), ]

# Maxlags=4 to capture annual effects.
# type="const" includes a constant/intercept in the VAR model.
cat("\n--- Step 2: Optimal VAR Lag Selection (up to maxlag=4) ---\n")
lag_select <- VARselect(ts_data_levels, lag.max = 4, type = "const", exogen = exog_vars_aligned)

optimal_lag <- lag_select$selection["AIC(n)"]
cat(paste0("Optimal VAR Lag (p) based on AIC: ", optimal_lag, "\n"))
k_ar_diff <- optimal_lag - 1 # Lag for Johansen test (p-1)
cat(paste0("Lag for Johansen Test (p-1): ", k_ar_diff, "\n"))


# ----------------------------------------------------------------------------------
# --- Step 3: Perform the Cointegration Test (Johansen) ---
# Use the optimal lag found (p=4, k_ar_diff=3 in our previous Python run)

# model="const" assumes a constant in the cointegrating relation (det_order=0 in Python)
# HINT: The urca package for Johansen does not directly accept the 'exogen' argument
# for deterministic dummies. The best practice is to include them in the subsequent 
# VECM estimation, which is what we did implicitly in the Python analysis.
cat("\n--- Step 3: Johansen Cointegration Test (using optimal k_ar_diff) ---\n")

# Use 'cajorls' model (case 3) which assumes a constant in the cointegrating vector
# and no deterministic trend in the levels data.
# Note: The test uses 'k' for VAR lag, which is 'p' in standard notation.
# We use p=4 as found by AIC.
johansen_test <- ca.jo(ts_data_levels, type = "trace", K = optimal_lag, ecdet = "const") 
summary(johansen_test)

# Interpretation:
r0_trace_stat <- johansen_test@teststat[1]
r0_crit_val <- johansen_test@cval[1, 2] # 5% critical value
final_rank <- 0

if (r0_trace_stat > r0_crit_val) {
  # If we reject H0: r=0 (i.e., Trace Stat > Critical Value), then rank is at least 1
  final_rank <- 1
}

# ----------------------------------------------------------------------------------
# --- Final Model Decision ---
cat("\n--- FINAL MODEL DECISION ---\n")
if (final_rank == 1) {
  cat(paste0("Cointegration Rank is r=1.\n"))
  cat("Decision: Use VECM (Vector Error Correction Model) and include dummies as 'exogen'.\n")
} else {
  cat(paste0("Cointegration Rank is r=0.\n"))
  cat("Decision: Use VAR on Combined Differences (Delta_1 * Delta_4).\n")
  cat("The long-run relationship is absent. Proceed to Step 4.\n")
}


# ----------------------------------------------------------------------------------
# Load libraries (ensuring they are still loaded from previous run)
library(vars)
library(urca)

# --- VECM Estimation Parameters ---
# The co.vec parameter is the cointegration vector (i.e., the error correction term)
# The adjustment parameter is the loading matrix (how quickly log_p0 and log_p1 adjust)

# The ca.jo test result that produced r=1
jo_result <- trend_joh 

# Estimated Cointegration Rank
r_rank <- 1 

# --- Estimate the VECM using the ca.jo output ---
# cajorls estimates the VECM coefficients based on the Johansen test results
# The function automatically uses the correct lag (p=3) and includes the 
# deterministic terms (trend + seasonal dummies) used in the jo_result object.
cat("\n--- Step 4: Estimating VECM with Rank r=1 ---\n")
vecm_model <- cajorls(jo_result, r = r_rank)
print(summary(vecm_model))

# Extract the fitted VECM results for further diagnostics (Step 5)
# This model now includes the long-run error correction term (ECT)
# and all the short-run dynamics (lags 1-3).

# 1. Convert the VECM results back to a VAR object for standard diagnostic tests
# The 'vars' package has better diagnostic tools than 'urca'
vecm_to_var <- vec2var(jo_result, r = r_rank)

cat("\n--- Step 5: Diagnostic Tests on VECM Residuals ---\n")

# 2. Test for Residual Serial Correlation (Autocorrelation)
# Null Hypothesis (H0): No serial correlation up to lag 8.
# If p-value > 0.05, the model is adequate.
cat("\n--- 5a: Serial Correlation Test (Portmanteau) ---\n")
serial_test <- serial.test(vecm_to_var, lags.pt = 8, type = "PT.asymptotic")
print(serial_test)
cat("\nIf p-value < 0.05, the model is misspecified (try increasing lag p).\n")

# 3. Test for Residual Normality
# Null Hypothesis (H0): Residuals are normally distributed.
cat("\n--- 5b: Residual Normality Test (Jarque-Bera) ---\n")
norm_test <- normality.test(vecm_to_var, multivariate.only = TRUE)
print(norm_test)
cat("\nIf p-value < 0.05, the residuals are non-normal (common in GDP data).\n")

# 4. Test for Residual Heteroscedasticity (ARCH effects)
# Null Hypothesis (H0): No ARCH effects.
cat("\n--- 5c: Residual Heteroscedasticity Test ---\n")
arch_test <- arch.test(vecm_to_var, lags.single = 4, multivariate.only = TRUE)
print(arch_test)
cat("\nIf p-value < 0.05, the residuals exhibit ARCH effects (volatility clustering).\n")
