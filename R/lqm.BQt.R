#' \code{lqm.BQt} fits linear quantile regression model
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
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is 5000
#' @param quiet see rjags package
#'
#' @return A \code{BQt} object which is a list with the following elements:
#'    \describe{
#'   \item{\code{Coefficients}}{list of posterior median of each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95}
#'   \item{\code{data}}{data include in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters}
#'    }
#'
#' @author Antoine Barbieri
#'
#' @export
#'
#' @examples
#' #---- Use data
#' data(wave)
#'
#' #---- Fit regression model for the first quartile
#' BQt_025 <- lqm.BQt(formula = h110d~vent_vit_moy,
#'                    data = wave,
#'                    tau = 0.25)
#'
#' #---- Get the estimated coefficients
#' BQt_025$coefficients
#'
#' #---- Summary of output
#' summary.BQt(BQt_025)
#'
lqm.BQt <- function(formula,
                    data,
                    tau = 0.5,
                    n.chains = 1,
                    n.iter = 10000,
                    n.burnin = 5000,
                    n.thin = 5,
                    n.adapt = 5000,
                    quiet = FALSE){

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
  parms_to_save <- c("beta", "sigma")

  #---- use JAGS sampler
  # if (!require("rjags"))
  #   stop("'rjags' is required.\n")
  JMjags.model <- rjags::jags.model(file = "JagsModel.txt",
                                    data = jags.data,
                                    # inits = list(initial.values),
                                    n.chains = n.chains,
                                    n.adapt = n.adapt,
                                    quiet = quiet)
  update(JMjags.model, n.burnin)
  fit <- rjags::coda.samples(JMjags.model,
                             variable.names = parms_to_save,
                             n.iter = n.iter - n.burnin,
                             thin = n.thin)
  file.remove(file.path(working.directory, "JagsModel.txt"))

  #---- MCMClist management
  Bs <- do.call(rbind, fit)
  sims.list <- vector("list", length(parms_to_save))
  names(sims.list) <- parms_to_save
  for (p in seq_along(parms_to_save)) {
    ii <- grep(paste("^", parms_to_save[p], sep = ""), colnames(Bs))
    sims.list[[p]] <- Bs[, ii]
  }
  sims.list <- list(beta = sims.list$beta, sigma = sims.list$sigma)
  colnames(sims.list$beta) <- colnames(X)

  #---- output
  out <- list(data = data)
  out$control <- list(fit = fit,
                      sims.list = sims.list,
                      formula = formula,
                      tau =tau,
                      call_function = "lqm.BQt")

  out$CIs <- lapply(sims.list, function(x) if (is.matrix(x))
    apply(x, 2, quantile, probs = c(0.025, 0.975))
    else quantile(x, probs =  c(0.025, 0.975)))
  out$coefficients <- lapply(sims.list, function(x) if (is.matrix(x))
    apply(x, 2, median)
    else median(x))
  class(out) <- "BQt"
  out

  }
