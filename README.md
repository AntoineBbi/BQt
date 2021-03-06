# BQt package for R

BQt is a R-package dealing the quantile regression in Bayesian framework. Based on the asymmetric Laplace distribution, it allows to estimate joint models for longitudinal and time-to-event data, linear mixed effects models and simple linear models.

As the JMcuR R-package, BQt is built from functions and raw codes of the JMbayes package version 0.4-1 implemented by Dimitris Rizopoulos. Another reference used to develop this package is:

Yang, M., Luo, S., & DeSantis, S. (2019). Bayesian quantile regression joint models: Inference and dynamic predictions. Statistical Methods in Medical Research, 28(8), 2524–2537. https://doi.org/10.1177/0962280218784757

To try the current development version from github, use:

```{r} 
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")}
devtools::install_github("AntoineBbi/BQt")
 ```
**Warning:** BQt package requires JAGS software (http://mcmc-jags.sourceforge.net/). 
