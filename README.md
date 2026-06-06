# Dynamic Analysis of Agriculture and GDP Growth in Malaysia

## Overview

This project investigates the dynamic relationship between Malaysia's agricultural sector and aggregate economic growth using multivariate time-series econometric techniques.

Quarterly data from 2015Q1 to 2025Q2 are analysed using Vector Autoregressive (VAR) models to evaluate whether agriculture contributes to GDP growth and how economic shocks propagate through the system over time.

The analysis combines forecasting, causality testing, impulse response analysis, and variance decomposition to examine both predictive and structural relationships.

---

## Research Questions

1. Does agricultural growth contribute to GDP growth in Malaysia?

2. Is there evidence of long-run equilibrium between agriculture and GDP?

3. Can agricultural activity predict future economic performance?

4. How do agricultural shocks affect GDP over time?

5. Which VAR specification provides superior forecasting performance?

---

## Dataset

| Feature | Description |
|----------|-------------|
| Country | Malaysia |
| Frequency | Quarterly |
| Period | 2015Q1 – 2025Q2 |
| Variables | GDP, Agricultural Output |
| Source | DOSM |

---

## Methodology

### Time Series Preprocessing

- Log Transformation
- Seasonal Adjustment
- Differencing

### Stationarity Analysis

- ADF Test
- KPSS Test
- HEGY Seasonal Unit Root Test
- Canova-Hansen Test

### Long-Run Analysis

- Johansen Cointegration Test
- Error Correction Testing

### Dynamic Modelling

- VAR Models
- Lag Selection
- Granger Causality

### Dynamic Interpretation

- Impulse Response Functions (IRF)
- Forecast Error Variance Decomposition (FEVD)

### Forecasting

- RMSE
- MAE
- MAPE
- Theil's U

---

## Key Findings

### Agriculture Granger-Causes GDP

Evidence suggests an agriculture-led growth channel where agricultural activity helps predict future GDP movements.

### Weak Cointegration Evidence

Long-run equilibrium relationships were statistically fragile, supporting the use of VAR rather than VECM.

### Positive Agriculture Shocks Support GDP

Impulse response analysis showed that positive agricultural shocks generate temporary increases in GDP growth.

### Agriculture Explains Part of GDP Variability

Forecast error variance decomposition indicates agriculture contributes approximately 17–23% of GDP variation over longer horizons.

### Forecast Performance Exceeds Naive Benchmarks

VAR models provided useful forecasting improvements for agricultural output and selected macroeconomic variables.

---

## Technical Skills Demonstrated

### Econometrics

- Time Series Econometrics
- VAR Modelling
- Cointegration Analysis
- Granger Causality
- Forecast Evaluation

### Statistics

- Hypothesis Testing
- Model Diagnostics
- Time Series Forecasting

### Programming

- R
- Forecasting Libraries
- Data Visualisation

---

## Technologies Used

- R
- vars
- urca
- forecast
- ggplot2

---
