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
      cat(paste("        - Correlation structure:  'middle' (autoregressive structure based on quantiles'order distance) \n"))
    if(object$control$corr_structure=="none")
      cat(paste("        - Correlation structure of b: only through covariance matrix of random effects  \n"))
    if(object$control$RE_ind)
      cat(paste("        - Correlation structure: quantile-specific covariance matrix of random effects  \n"))
    if(!object$control$RE_ind)
      cat(paste("        - Correlation structure of b: common covariance matrix of random effects  \n"))
  }
  cat(paste("     - Quantile order(s): ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
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
    if(object$control$call_function %in% c("lqm.BQt","lqmm.BQt","qrjm.BQt")){
      beta_estim <- cbind("Value" = coefs$beta,
                           "2.5%" = CIs$beta[1, ],
                           "97.5%" = CIs$beta[2, ],
                           "Rhat" = Rhat$beta)
      sigma_estim <- cbind("Value" = coefs$sigma,
                           "2.5%" = CIs$sigma[1],
                           "97.5%" = CIs$sigma[2],
                           "Rhat" = Rhat$sigma)
      rownames(sigma_estim) <- "sigma"
    }else{
      beta_estim <- cbind("Value" = as.vector(t(coefs$beta)),
                          "2.5%" = CIs$beta[, 1],
                          "97.5%" = CIs$beta[, 2],
                          "Rhat" = as.vector(t(Rhat$beta)))
      sigma_estim <- cbind("Value" = coefs$sigma,
                           "2.5%" = CIs$sigma[1],
                           "97.5%" = CIs$sigma[2],
                           "Rhat" = Rhat$sigma)
      rownames(sigma_estim) <- names(coefs$sigma)
    }
    cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
    prmatrix(beta_estim, na.print = "")
    cat("\n")
    cat("#-- Estimation of 'sigma' parameter associated with asymmetric Laplace distribution: \n")
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
    # for multiple quantile
    if(object$control$call_function %in% c("mlqmm.BQt","mqrjm.BQt")){
      if(object$control$RE_ind){
        tmp <- round(rbind(cbind(object$mean$covariance.b1,
                           matrix(0, nrow=nrow(object$mean$covariance.b1),
                                  ncol=ncol(object$mean$covariance.b2)+ncol(object$mean$covariance.b3))),
                     cbind(matrix(0,
                                  nrow=nrow(object$mean$covariance.b2),
                                  ncol=ncol(object$mean$covariance.b1)),
                           object$mean$covariance.b2,
                           matrix(0,
                                  nrow=nrow(object$mean$covariance.b2),
                                  ncol=ncol(object$mean$covariance.b3))),
                     cbind(matrix(0,
                                  nrow=nrow(object$mean$covariance.b3),
                                  ncol=ncol(object$mean$covariance.b1)+ncol(object$mean$covariance.b2)),
                           object$mean$covariance.b3)),
                     3)
        colnames(tmp) <- rownames(tmp) <- paste(rep(paste("tau", as.character(object$control$tau*100), sep = ""),
                                                    each = nrow(object$mean$covariance.b1)),
                                                rep(colnames(object$mean$covariance.b1), 3),
                                                sep = ".")
        prmatrix(tmp, na.print = "")

      }
      if(!object$control$RE_ind)
        prmatrix(object$mean$covariance.b, na.print = "")
      if(object$control$corr_structure!="none"){
        if(length(object$control$tau)==3 & object$control$corr_structure=="free"){
          rho_estim <- cbind("Value" = c(coefs$rho12, coefs$rho13, coefs$rho23),
                              "2.5%" = c(CIs$rho12[1], CIs$rho13[1], CIs$rho23[1]),
                              "97.5%" = c(CIs$rho12[2], CIs$rho13[2], CIs$rho23[2]),
                              "Rhat" = c(Rhat$rho12, Rhat$rho13, Rhat$rho23))
          rownames(rho_estim) <- c("corr.Q1.Q2", "corr.Q1.Q3", "corr.Q2.Q3")
        }else{
          rho_estim <- cbind("Value" = as.vector(coefs$rho),
                             "2.5%" = as.vector(CIs$rho[1]),
                             "97.5%" = as.vector(CIs$rho[2]),
                             "Rhat" = as.vector(Rhat$rho))
          rownames(rho_estim) <- "rho"
        }
        cat("\n")
        cat("#-- 'rho' parameter(s) associated with the correlation strucutre : \n")
        prmatrix(rho_estim, na.print = "")
      }
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
