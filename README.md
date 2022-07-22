# financial-risk-measurement
This repo contains the source code where I implement the topics of the Risk Econometrics exam that I studied during my academic experience at the University of Bologna.

The textbook from which I studied the theory is *"[Elements of Financial Risk Management](https://www.amazon.it/Elements-Financial-Management-Peter-Christoffersen/dp/0128102357)"*, and the implementation was conducted in [R](https://cran.r-project.org/), using the `rugarch`,` rmgarch` and `quarks` packages.

## Goal
The goal of the project is to measure the portfolio risk using [Value at Risk](https://en.wikipedia.org/wiki/Value_at_risk) and the [Expected Shortfall](https://en.wikipedia.org/wiki/Expected_shortfall), comparing the different methods for their calculation and modeling the conditional volatility of log returns.

The methods used are applied separately considering both a univariate and a multivariate portfolio.

## Univariate analysis
In the first part of the project, [`univariate-financial-risk-measurement`](https://github.com/gianlucaciaccio/financial-risk-measurement/tree/main/univariate-financial-risk-measurement), is considered a portfolio with a single asset, the SP500 from January 1st 2001 to December 31st 2010. 


The methods used are:

#### **Case with unconditional (static) volatility over time:**

- ***Non-parametric methods*** (no assumptions on the distribution of returns)

  - *Historical Simulation (HS)*
  - *Age-Weighted Historical Simulation (WHS)*

- ***Parametric methods*** (hypothesis on the distribution of returns)
  - *Normal Distribution*
  - *Stundent-t distribution*
  - *Cornish-Fisher expansion*

VaR and ES are calculated both on the entire portfolio and using rolling windows of 250 and 1000 trading days.

#### **Case with conditional (dynamic) volatility over time:**

Different methods for modeling conditional volatility are compared: 

- The RiskMetrics *EWMA model*
- Different specifications of *GARCH* (*Standard GARCH*, *Integrated-GARH*, *GJR-GARCH*, *Nonlinear-GARCH* and *Exponential-GARCH*).

For each of these models, different distributional assumptions of returns have been evaluated:

- *Normal distribution*
- *Stundent-t distribution*
- Empirical distribution using the *Cornish-Fisher expansion* method
- Empirical distribution using the *Filtered Historical Simulation (FHS)* method
- Empirical distribution using the *Extreme Value Theory (EVT)* method (Hill estimator)

The best specification was chosen using *backtesting*.


## Multivariate analysis
In the second part of the project, [`multivariate-financial-risk-measurement`](https://github.com/gianlucaciaccio/financial-risk-measurement/tree/main/multivariate-financial-risk-measurement), is considered a portfolio with two assets, one equity (SP500) and one bond (U.S. 10 Year Treasury Note). 

By maintaining a similar scheme of the univariate analysis, the focus is only on modeling portfolio volatility and covariance/correlation between assets over time using *Dynamic Conditional Correlation* models (*DCC-EWMA* and *DCC-GARCH*).
