#' Summary of a \code{BQt} objects
#'
#' The function provides a summary of linear quantile regression model (LQRM) estimations. First, LQRM's details are given.
#' Next, the estimation parameters (posterior median) and bound of the credible interval at 95% are given.
#'
#' @param object an object inheriting from class 'JMcuR'
#' @param ... further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @return Returns NULL.
#'
#' @author Antoine Barbieri
#'
#' @seealso \code{\link{lqm.BQt}}
#'
summary.BQt <- function (object, ...)
  {
  if (!inherits(object, "BQt"))
    stop("Use only with 'BQt' objects.\n")

  #---- Global details
  cat("#-- Statistical model:", "\n")
  cat(paste("     - Quantile regression from '", object$control$call_function), "' function \n")
  if(object$control$call_function=="qrjm.BQt"){
    if(object$control$param=="sharedRE")
      cat(paste("        - Association structure in joint model: Shared random effect(s) \n"))
    if(object$control$param=="value")
      cat(paste("        - Association structure in joint model: Current quantile of longitudinal process \n"))
    if(object$control$survMod=="weibull")
      cat(paste("        - Baseline risk function in survival model: Weibull distribution \n"))
  }
  cat(paste("     - Quantile order: ", object$control$tau), "\n")
  cat(paste("     - Number of observations: ", nrow(object$data)), "\n")
  if(object$control$call_function %in% c("lqmm.BQt","qrjm.BQt")){
    cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  }
  # if(object$control$call_function=="qrjm.BQt"){
  #   cat(paste("     - Number of observed events (e.g. subject): ", sum(object$data$event), "\n"))
  # }
  cat("\n")

    #---- Parameter Estimations
    coefs <- object$coefficients
    CIs <- object$CIs
    # beta regression parameters
    param_estim <- cbind("Value" = coefs$beta,
                       "2.5%" = CIs$beta[1, ],
                       "97.5%" = CIs$beta[2, ])
    cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
    prmatrix(param_estim, na.print = "")
    cat("\n")
    cat("#-- Estimation of sigma parameter associated with asymmetric Laplace distribution: \n")
    sigma_estim <- cbind("Value" = coefs$sigma,
                         "2.5%" = CIs$sigma[1],
                         "97.5%" = CIs$sigma[2])
    rownames(sigma_estim) <- "sigma"
    prmatrix(sigma_estim,
             na.print = "")

    # Random effects for "mixed regression model
    if(object$control$call_function %in% c("lqmm.BQt","qrjm.BQt")){
      cat("\n")
      cat("#-- (Co)variance matrix of the random-effect(s): \n")
      prmatrix(object$coefficients$Sigma2, na.print = "")
    }

    # alpha regression parameters
    if(object$control$call_function=="qrjm.BQt"){
      cat("\n")
      param_estim <- cbind("Value" = coefs$alpha,
                           "2.5%" = CIs$alpha[1, ],
                           "97.5%" = CIs$alpha[2, ])
      cat("#-- Estimation of survival parameters and their credible interval bounds: \n")
      prmatrix(param_estim, na.print = "")
    }

    # survival structure parameter
    if(object$control$call_function=="qrjm.BQt" && object$control$survMod=="weibull"){
      cat("\n")
      param_estim <- cbind("Value" = coefs$shape,
                           "2.5%"  = CIs$shape[1],
                           "97.5%" = CIs$shape[2])
      cat("#-- Estimation of shape parameters and its credible interval bounds for Weibull model: \n")
      prmatrix(param_estim, na.print = "")
    }

  }
