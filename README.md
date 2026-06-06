# Agriculture and GDP Dynamics in Malaysia

This project investigates the dynamic relationship between Malaysia's agriculture sector and total GDP using multivariate time-series econometric methods.

Using quarterly data from 2015Q1 to 2025Q2, the study evaluates whether agricultural activity contributes to economic growth and how shocks propagate through the economy over time.

The analysis follows a complete econometric workflow including stationarity testing, seasonal unit root analysis, cointegration testing, VAR model estimation, diagnostic checking, forecasting, Granger causality analysis, impulse response functions (IRFs), and forecast error variance decomposition (FEVD).

The objective is to quantify both short-run dynamics and predictive relationships between agriculture and aggregate economic performance.

## Methodology

The modelling pipeline consists of:

1. Logarithmic transformation
2. Seasonal and trend decomposition
3. Unit root testing
   - ADF
   - KPSS
   - HEGY seasonal unit root test
   - Canova-Hansen seasonal stability test

4. Cointegration analysis
   - Johansen Trace Test
   - Error Correction Term stationarity testing

5. VAR model specification
   - Lag order selection using AIC, BIC, HQ and FPE
   - VAR with exogenous variables

6. Model diagnostics
   - Portmanteau test
   - ARCH-LM test
   - Jarque-Bera test
   - Stability analysis

7. Forecasting
   - Out-of-sample evaluation
   - RMSE
   - MAE
   - MAPE
   - Theil's U

8. Dynamic analysis
   - Granger causality
   - Impulse Response Functions (IRFs)
   - Forecast Error Variance Decomposition (FEVD)
