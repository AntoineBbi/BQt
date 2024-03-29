#' \code{lqm} fits linear quantile regression model
#'
#' Function using JAGS to estimate the linear quantile regression model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formula formula for the quantile regression including response variable
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is NULL
#' @param save_jagsUI If TRUE (is FALSE by default), the output of jagsUI package is return by the function
#' @param parallel see jagsUI::jags() function
#'
#' @return A \code{BQt} object which is a list with the following elements:
#'    \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{Rhat}}{Gelman and Rubin diagnostic for all parameters}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{W}}{list including both posterior mean and posterior standard deviation of subject-specific random varaible W}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#'
#' @author Antoine Barbieri
#'
#' @import jagsUI
#' @export
#'
#' @examples
#' #---- Use data
#' data(wave)
#'
#' #---- Fit regression model for the first quartile
#' lqm_025 <- lqm(formula = h110d~vent_vit_moy,
#'                data = wave,
#'                n.iter = 1000,
#'                n.burnin = 500,
#'                tau = 0.25)
#'
#' #---- Get the posterior mean of parameters
#' lqm_025$mean
#'
#' #---- Summary of output
#' summary(lqm_025)
#'
lqm <- function(formula,
                data,
                tau = 0.5,
                n.chains = 3,
                n.iter = 10000,
                n.burnin = 5000,
                n.thin = 1,
                n.adapt = NULL,
                save_jagsUI = FALSE,
                parallel = FALSE){

  # data
  tmp <- data[c(all.vars(formula))]
  tmp <- unique(tmp)
  I <- nrow(tmp)
  # design matrices
  y <- data[all.vars(formula)][,1]
  mfX <- model.frame(formula, data = tmp)
  X <- model.matrix(formula, mfX)
  ncX <- ncol(X)
  # use lm function to initiated values
  lm_tmp <- stats::lm(formula,
                      data = tmp)
  # prior beta parameters
  priorMean.beta <- as.numeric(lm_tmp$coefficients)
  # priorTau.beta <- diag(1/c(1,smcure_out$beta_var)/100)
  # list of data jags
  jags.data <- list(y = y,
                    X = X,
                    tau = tau,
                    ncX = ncX,
                    I = I
                    # ,
                    # priorMean.beta = priorMean.beta,
                    # priorTau.beta = priorTau.beta,
                    # priorA.shape = 1/100,
                    # priorB.shape = 1/100
                    )

  #---- write jags model in txt from R function
  working.directory = getwd()
  write.model.jags(model = jags_lqm,
                   name_model = "jags_lqm",
                   intitled = file.path(working.directory,"JagsModel.txt"),
                   Data = jags.data)

  # jags argument
  parms_to_save <- c("beta", "sigma", "va1")

  #---- use JAGS sampler
  # if (!require("rjags"))
  #   stop("'rjags' is required.\n")

  # using jagsUI
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = "JagsModel.txt",
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          DIC = F)

  file.remove(file.path(working.directory, "JagsModel.txt"))

  #---- output
  out <- list(data = data)
  out$control <- list(formula = formula,
                      tau = tau,
                      call_function = "lqm",
                      n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin)

  #- other outputs

  # sims.list output
  out$sims.list <- out_jags$sims.list
  out$sims.list$va1 <- NULL

  # Random variable W esponentially distributed
  out$W <- list(postMeans = out_jags$mean$va1,
                postSd = out_jags$sd$va1)

  # median : Posterior median of parameters (if mean, you can use mean instead of q50)
  out$median <- out_jags$q50
  out$median$va1 <- NULL

  # mean : Posterior mean of parameters (if mean, you can use mean instead of q50)
  out$mean <- out_jags$mean
  out$mean$va1 <- NULL

  # modes of parameters
  out$modes <- lapply(out$sims.list, function(x) {
    m <- function(x) {
      d <- density(x, bw = "nrd", adjust = 3, n = 1000)
      d$x[which.max(d$y)]
    }
    if (is.matrix(x))
      as.array(apply(x, 2, m))
    else{
      if(is.array(x))
        apply(x, c(2,3), m)
      else m(x)
    }
  })

  # standard error of parameters
  out$StErr <- lapply(out$sims.list, function(x) {
    f <- function(x) {
      acf.x <- drop(acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
      acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
      ess <- length(x)/(1 + 2 * sum(acf.x))
      sqrt(var(x)/ess)
    }
    if (is.matrix(x))
      as.array(apply(x, 2, f))
    else{
      if(is.array(x))
        apply(x, c(2,3), f)
      else f(x)
    }
  })

  # standard deviation of parameters
  out$StDev <- out_jags$sd
  out$StDev$va1 <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  out$Rhat$va1 <- NULL

  # names
  names(out$mean$beta) <-
    names(out$median$beta) <-
    names(out$Rhat$beta) <-
    names(out$modes$beta) <-
    names(out$StErr$beta) <-
    names(out$StDev$beta) <- colnames(X)

  # credible intervalles
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- colnames(X)
  colnames(out$CIs$beta) <- c("2.5%", "97.5%")

  out$CIs$sigma <- c(out_jags$q2.5$sigma,
                     out_jags$q97.5$sigma)
  names(out$CIs$sigma) <- c("2.5%", "97.5%")

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  # End of function
  class(out) <- "BQt"
  out

  }
