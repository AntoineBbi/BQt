#' Summary of a \code{BQt} objects
#'
#' The function provides a summary of BQt object (lqm, lqmm, qrjm, mlqmm).
#' Therfore, the estimation parameters (posterior mean) and bounds of the credible interval at 95% and Gelman & Rubin diagnostic are given.
#'
#' @param object an object inheriting from class 'BQt'
#' @param ... further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @return Returns NULL.
#'
#' @author Antoine Barbieri
#'
#' @seealso \code{\link{lqm.BQt}}, \code{\link{lqmm.BQt}}, \code{\link{mlqmm.BQt}}, \code{\link{qrjm.BQt}}
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
  if(object$control$call_function=="mlqmm.BQt"){
    if(object$control$corr_structure=="free")
      cat(paste("        - Correlation structure:  'Free' \n"))
    if(object$control$corr_structure=="middle")
      cat(paste("        - Correlation structure:  'middle' autoregressive structure based on quantiles'order distance \n"))
    if(object$control$corr_structure=="none")
      cat(paste("        - Correlation structure: only through covariance matrix of random effects  \n"))
  }
  cat(paste("     - Quantile order(s): ", object$control$tau), "\n")
  cat(paste("     - Number of observations: ", nrow(object$data)), "\n")
  if(object$control$call_function %in% c("lqmm.BQt","qrjm.BQt", "mlqmm.BQt")){
    cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  }
  if(object$control$call_function=="qrjm.BQt"){
    cat(paste("     - Number of observed events: ", sum(object$control$event), "\n"))
  }
  cat("\n")

    #---- Parameter Estimations
    coefs <- object$mean
    CIs <- object$CIs
    Rhat <- object$Rhat
    # beta regression parameters
    param_estim <- cbind("Value" = coefs$beta,
                       "2.5%" = CIs$beta[1, ],
                       "97.5%" = CIs$beta[2, ],
                       "Rhat" = Rhat$beta)
    cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
    prmatrix(param_estim, na.print = "")
    cat("\n")
    cat("#-- Estimation of sigma parameter associated with asymmetric Laplace distribution: \n")
    sigma_estim <- cbind("Value" = coefs$sigma,
                         "2.5%" = CIs$sigma[1],
                         "97.5%" = CIs$sigma[2],
                         "Rhat" = Rhat$sigma)
    rownames(sigma_estim) <- "sigma"
    prmatrix(sigma_estim, na.print = "")

    # Random effects for "mixed regression model
    if(object$control$call_function %in% c("lqmm.BQt","qrjm.BQt")){
      cat("\n")
      cat("#-- (Co)variance matrix of the random-effect(s): \n")
      if(object$control$RE_ind){
        tmp <- diag(object$mean$covariance.b)
        colnames(tmp) <- rownames(tmp) <- names(object$mean$covariance.b)
        prmatrix(tmp, na.print = "")

      }
      if(!object$control$RE_ind)
        prmatrix(object$mean$covariance.b, na.print = "")
    }

    # survival parameters
    # alpha regression parameters
    if(object$control$call_function=="qrjm.BQt"){
      # survival structure parameter
      if(object$control$survMod=="weibull"){
        param_estim1 <- cbind("Value" = object$mean$shape,
                              "2.5%"  = object$CIs$shape[1],
                              "97.5%" = object$CIs$shape[2],
                              "Rhat" = object$Rhat$shape)
        rownames(param_estim1) <- "shape"
      }
      # alpha parameters
      param_estim2 <- cbind("Value" = object$mean$alpha,
                            "2.5%" = object$CIs$alpha[, 1],
                            "97.5%" = object$CIs$alpha[, 2],
                            "Rhat" = object$Rhat$alpha)
      # association parameters
      if(length(object$mean$alpha.assoc)==1)
        param_estim3 <- cbind("Value" = object$mean$alpha.assoc,
                              "2.5%" = object$CIs$alpha.assoc[1],
                              "97.5%" = object$CIs$alpha.assoc[2],
                              "Rhat" = object$Rhat$alpha.assoc)
      else
        param_estim3 <- cbind("Value" = object$mean$alpha.assoc,
                              "2.5%" = object$CIs$alpha.assoc[, 1],
                              "97.5%" = object$CIs$alpha.assoc[, 2],
                              "Rhat" = object$Rhat$alpha.assoc)
      rownames(param_estim3) <- rep("alpha.assoc", length(object$mean$alpha.assoc))
      # Print
      cat("\n")
      cat("#-- Estimation of survival models: \n")
      prmatrix(rbind(param_estim1, param_estim2, param_estim3), na.print = "")
    }

  }
