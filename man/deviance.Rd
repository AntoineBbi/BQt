% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deviance.R
\name{deviance}
\alias{deviance}
\title{\code{deviance} returns the deviance based on the conditional likelihood associated with the survival part.}
\usage{
deviance(object, M = 1000, conditional = "survival", verbose = TRUE)
}
\arguments{
\item{object}{an object inheriting from class 'BQt'.}

\item{M}{an integer indicating the number of draws used for the approximation of the integral with respect to random effects, M=1000 by default.}

\item{conditional}{is "survival" by default because only this one is implemented until now.}

\item{verbose}{A logical indicating if information about method's progress (included progress bars for each step) must be printed (default to TRUE). Adds a small extra overload.}
}
\value{
An object which is a list with the following elements:
   \describe{
  \item{\code{deviance}}{Numerical object returning the deviance}
  \item{\code{likelihood}}{(Conditional) likelihood}
  \item{\code{sims.list}}{list of individual quantities like likelihood, draws of random effects, hazard and survival functions}
  \item{\code{control}}{list of arguments giving details about the deviance}
 }
}
\description{
\code{deviance} returns the deviance based on the conditional likelihood associated with the survival part.
}
\examples{

\dontrun{
data("aids", package = "joineR")

#---- Fit quantile regression joint model for the first quartile
qrjm_5 <- qrjm.BQt(formFixed = CD4 ~ obstime,
                   formRandom = ~ obstime,
                   formGroup = ~ id,
                   formSurv = Surv(time, death) ~ drug + gender + prevOI + AZT,
                   survMod = "weibull",
                   n.iter = 1000,
                   n.burnin = 500,
                   n.thin = 1,
                   n.adapt = 500,
                   param = "value",
                   timeVar= "obstime", save_va = TRUE, parallel = TRUE,
                   data = aids,
                   tau = 0.5)

deviance(qrjm_5, M=200)
}

}
\author{
Antoine Barbieri and Baptiste Courrèges
}
