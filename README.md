# financial-risk-measurement
This repo contains the source code of my project carried out during my academic experience at the University of Bologna.

The goal of the project is to measure the portfolio risk using Value at Risk and the Expected Shortfall, comparing the different methods for their calculation and focusing on modeling the conditional volatility of log returns.

In particular, in the first part of the project, **"Univariate financial risk measurement"** is considered a portfolio with a single asset (the SP 500 index on a time series from January 1st 2001 to December 31st 2010). The calculation of VaR and ES is first done without considering the conditional nature of volatility over time using the non-parametric method of Historical Simulation and its Age-Weighted version, using floating windows of 250 and 1000 trading days. Alongside this approach, the two risk measures are adapted assuming a normal and Stundent-t distribution of the log returns and using the Cornish-Fisher approximation. 

Subsequently, volatility is modeled through the approach proposed by RiskMetrics and the different GARCH specifications (Integrated-GARH, GJR-GARCH, Nonlinear-GARCH and Exponential-GARCH). For these models, the best specification is chosen through VaR backtesting.

In the second part of the project, **"Multivariate Financial Risk Measurement"**, we consider a portfolio with two assets, one equity (SP 500) and one bond (U.S. 10 Year Treasury Note). By maintaining the same scheme of the univariate analysis, in this phase the portfolio risk is measured considering also the dynamic covariance and correlation of the two assets.
