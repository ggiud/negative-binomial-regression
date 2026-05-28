readme = """# COVID-19 Statistical Analysis & Count Models (ECDC Data)

## Overview

This project analyzes COVID-19 weekly death counts in Europe (2020–2023) using official ECDC data.

The goal is to:
- Understand the statistical structure of pandemic data
- Evaluate count regression models
- Handle overdispersion and zero inflation
- Compare predictive performance across models

---

## Dataset

Source: European Centre for Disease Prevention and Control (ECDC)

- Time range: 2020–2023 (weekly data)
- Geography: European countries
- Observations: ~12,600
- Key variables: cases, deaths, cumulative counts, incidence rates

Target variable:
weekly_count_deaths

---

## Preprocessing

- Removed non-informative variables (country_code, continent, source, note)
- Split data by indicator (cases / deaths)
- Removed missing values (~5.5%)
- Final dataset: 9 variables

---

## Exploratory Data Analysis

- Strong positive correlations between variables
- Highest correlation: cumulative cases vs deaths
- Presence of outliers
- Heavy skewness in distributions
- Target variable does NOT follow a Poisson distribution

---

## Models

### Poisson Regression
- Poor fit
- Severe overdispersion
- Underestimated standard errors

### Negative Binomial Regression
- Handles overdispersion well
- Best classical GLM baseline

### Alternative Models
- Quasi-Poisson
- Zero-Inflated Poisson (ZIP)
- Zero-Inflated Negative Binomial (ZINB)
- Hurdle Model

Best performance:
- ZINB (prediction)
- Hurdle model (explanatory power)

---

## Model Comparison

Metrics used:
- MSE
- RMSE
- MAE
- R²
- AIC / BIC

Key findings:
- Poisson performs worst
- Negative Binomial is strong baseline
- ZINB performs best overall
- Hurdle model has highest R²

---

## Conclusion

COVID-19 death counts exhibit:
- Overdispersion
- Zero inflation
- Poor Poisson fit

Best models:
- ZINB (predictive performance)
- Negative Binomial (robust baseline)

"""

