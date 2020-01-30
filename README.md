# BQt package for R

BQt is a R-package dealing the quantile regression in Bayesian framework. Based on the asymmetric Laplace distribution, it allows to estimate joint models for longitudinal and time-to-event data, linear mixed effects models and simple linear models.

To try the current development version from github, use:

```{r} 
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")}
devtools::install_github("AntoineBbi/BQt")
 ```
**Warning:** BQt package requires JAGS software (http://mcmc-jags.sourceforge.net/). 
