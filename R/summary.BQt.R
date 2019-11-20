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
  cat("#-- Statistical Model:", "\n")
  cat(paste("     - Quantile regression from '", object$call_function), "' function \n")
  cat(paste("     - Quantile order:", object$tau), "\n")
  cat(paste("     - Number of observations:", nrow(object$data)), "\n")
  cat("\n")

    #---- Parameter Estimations
    coefs <- object$coefficients
    CIs <- object$CIs
    # beta regression parameters
    param_estim <- cbind("Value" = coefs$beta,
                       "2.5%" = CIs$beta[1, ],
                       "97.5%" = CIs$beta[2, ])
    cat("#-- Estimation of regression parameters: \n")
    prmatrix(param_estim, na.print = "")
    cat("\n")
    cat("#-- Estimation of sigma parameter associated with asymmetric Laplace distribution: \n")
    sigma_estim <- cbind("Value" = coefs$sigma,
                         "2.5%" = CIs$sigma[1],
                         "97.5%" = CIs$sigma[2])
    rownames(sigma_estim) <- "sigma"
    prmatrix(sigma_estim,
             na.print = "")

  }
