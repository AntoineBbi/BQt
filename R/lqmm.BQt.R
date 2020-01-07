#' \code{lqmm.BQt} fits linear quantile mixed model
#'
#' Function using JAGS to estimate the linear quantile mixed model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param RE_ind Boolean denoting if the random effects are assumed independent ; default is \code{FALSE}
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is 5000
#' @param quiet see rjags package
#' @param precision variance by default for vague prior distribution
#'
#' @return A \code{BQt} object is a list with the following elements:
#'  \describe{
#'   \item{\code{coefficients}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters excepted random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{postMeanq}}{data including posterior mean of subject-specific random effects}
#'   \item{\code{postVars}}{list of subject-specific random effect covariance matrix}
#'  }
#'
#
#' @author Antoine Barbieri
#'
#' @import rjags coda lqmm
#'
#' @references Marco Geraci and Matteo Bottai (2014).
#' \emph{Linear quantile mixed models}.
#' Statistics and Computing, 24(3):461-479. doi: 10.1007/s11222-013-9381-9.
#'
#' @export
#'
#' @examples
#'
#' #---- Orthodont data from lqmm package
#' data("Orthodont", package = "lqmm")
#'
#' #---- Fit regression model for the first quartile
#' BQt_025 <- lqmm.BQt(formFixed = distance ~ age,
#'                     formRandom = ~ age,
#'                     formGroup = ~ Subject,
#'                     data = Orthodont,
#'                     tau = 0.25)
#'
#' #---- Get the estimated coefficients
#' BQt_025$coefficients
#'
#' #---- Summary of output
#' summary.BQt(BQt_025)
#'
lqmm.BQt <- function(formFixed,
                     formRandom,
                     formGroup,
                     data,
                     tau,
                     RE_ind = FALSE,
                     n.chains = 1,
                     n.iter = 10000,
                     n.burnin = 5000,
                     n.thin = 5,
                     n.adapt = 5000,
                     quiet = FALSE,
                     precision = 10){

  # To do: improve chain convergence.

  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  # use lqmm function to initiated values
  cat("Initiation of parameter values using lqmm package. \n")
  tmp_model <- lqmm::lqmm(fixed = formFixed,
                          random = formRandom,
                          group = id,
                          tau = tau,
                          data = data_long)
  # prior beta parameters
  priorMean.beta <- coef(tmp_model)
  priorTau.beta <- diag(rep(1/10,length(priorMean.beta)))

  bis <- as.matrix(ranef(tmp_model))
  bis[abs(bis)<.0001] <- 0
  initial.values <- list(b = bis,
                         beta = priorMean.beta,
                         sigma = tmp_model$scale)

  # list of data jags
  jags.data <- list(y = y,
                    X = X,
                    U = U,
                    tau = tau,
                    ncX = ncol(X),
                    ncU = ncol(U),
                    I = I,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorA.sigma = 1/precision,
                    priorB.sigma = 1/precision
                    )

  if(jags.data$ncU==1)
    RE_ind <- TRUE
  if(RE_ind){
    jags.data <- c(jags.data,
                   list(priorA.Sigma2 = 1/precision,
                        priorB.Sigma2 = 1/precision)
                   )
    initial.values$prec.Sigma2 <- 1/VarCorr(tmp_model)
  }else{
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(1/precision, ncol(U))),
                        priorK.Sigma2 = ncol(U),
                        mu0 = rep(0, ncol(U))
                        )
                   )
    initial.values$prec.Sigma2 <- diag(1/VarCorr(tmp_model))
  }
  model <- switch(paste(RE_ind, sep = "/"),
                  # wishart distribution when more than one RE is considered
                  `FALSE` = jags_lqmm,
                  # inverse gamma distribution when only one RE is considered
                  `TRUE` = jags_lqmm_IG
                  )

  # parameters to save in the sampling step
  parms_to_save <- c("beta", "sigma", "b", "prec.Sigma2")

  #---- write jags model in txt from R function
  working.directory = getwd()
  write.model.jags(model = model,
                   name_model = "jags_lqmm",
                   intitled = file.path(working.directory,"JagsModel.txt"),
                   Data = jags.data)

  #---- use JAGS sampler
  # if (!require("rjags"))
  #   stop("'rjags' is required.\n")
  JMjags.model <- rjags::jags.model(file = "JagsModel.txt",
                                    data = jags.data,
                                    inits = initial.values,
                                    n.chains = n.chains,
                                    n.adapt = n.adapt,
                                    quiet = quiet)
  update(JMjags.model, n.burnin)
  fit <- rjags::coda.samples(JMjags.model,
                             variable.names = parms_to_save,
                             n.iter = n.iter - n.burnin,
                             thin = n.thin)
  file.remove(file.path(working.directory, "JagsModel.txt"))

  #---- output building

  #-- MCMClist management

  #- arguments
  out <- list(data = data)
  out$control <- list(fit = fit,
                      formFixed,
                      formRandom,
                      formGroup,
                      tau = tau,
                      call_function = "lqmm.BQt",
                      I = I)

  #- other outputs
  Bs <- do.call(rbind, fit)
  sims.list <- vector("list", length(parms_to_save))
  names(sims.list) <- parms_to_save
  for (p in seq_along(parms_to_save)) {
    ii <- grep(paste("^", parms_to_save[p], sep = ""), colnames(Bs))
    sims.list[[p]] <- Bs[, ii]
  }
  sims.list <- list(beta = sims.list$beta,
                    sigma = sims.list$sigma,
                    b = sims.list$b,
                    prec.Sigma2 = sims.list$prec.Sigma2)
  colnames(sims.list$beta) <- colnames(X)

  #- random effect management
  ind.bs <- grep("b[", colnames(Bs), fixed = TRUE)
  sims.list$b <- sims.list$b[, ind.bs, drop = FALSE]
  ranef <- sims.list$b
  if(jags.data$ncU!=1){
    ord.col <- sapply(strsplit(colnames(ranef), "[", fixed = TRUE),"[", 2)
    ord.col <- sapply(strsplit(ord.col, ",", fixed = TRUE), "[",1)
    ord.col <- order(as.numeric(ord.col))
    ranef <- ranef[, ord.col, drop = FALSE]
  }
  postMeans <- matrix(colMeans(ranef), ncol = ncol(U), byrow = TRUE)
  dimnames(postMeans) <- list(levels(factor(id)), colnames(U))
  postVars <- vector("list", nrow(postMeans))
  ind.var <- matrix(seq_len(ncol(ranef)), ncol = ncol(U), byrow = TRUE)
  for (i in seq_along(postVars)){
    postVars[[i]] <- var(ranef[, ind.var[i, ]])
  }
  if(jags.data$ncU!=1){
    postVars[] <- lapply(postVars,
                         function(m) {
                           dimnames(m) <- list(colnames(postMeans), colnames(postMeans))
                           m
                         }
    )
  }
  names(postVars) <- rownames(postMeans)
  out$postMeans = postMeans
  out$postVars = postVars
  out$sims.list.b <- sims.list$b
  sims.list$b <- NULL

  #- draws of parameters
  out$sims.list <- sims.list
  # model parameters
  if(!RE_ind && jags.data$ncU>1){
    indSigma2 <- cbind(rep(1:jags.data$ncU, each = jags.data$ncU),
                       rep(1:jags.data$ncU, jags.data$ncU))
    Sigma2 <- t(sapply(seq_len(nrow(Bs)),
                       function(i) {
                         m <- matrix(0, jags.data$ncU, jags.data$ncU)
                         m[indSigma2] <- sims.list$prec.Sigma2[i, ]
                         d <- solve(m)
                         d[lower.tri(d, TRUE)]
                       }))
    out$sims.list$Sigma2 = Sigma2
    out$sims.list$prec.Sigma2 <- NULL
    tmp_mat <- matrix(1, ncol = jags.data$ncU, nrow = jags.data$ncU)
    colnames(out$sims.list$Sigma2) <- paste("Sigma2[", row(tmp_mat)[lower.tri(tmp_mat,TRUE)],
                                            ", ",
                                            col(tmp_mat)[lower.tri(tmp_mat, TRUE)],
                                            "]",
                                            sep = "")
  }else{ # if RE_ind
    out$sims.list$Sigma2 = 1/sims.list$prec.Sigma2
    out$sims.list$prec.Sigma2 <- NULL
    if(jags.data$ncU>1){
      colnames(out$sims.list$Sigma2) <- gsub("prec.", "", colnames(sims.list$Sigma2))
    }
  }
  # all parameters
  out$CIs <- lapply(out$sims.list, function(x) if (is.matrix(x))
    apply(x, 2, quantile, probs = c(0.025, 0.975))
    else quantile(x, probs =  c(0.025, 0.975)))
  out$coefficients <- lapply(out$sims.list, function(x) if (is.matrix(x))
    apply(x, 2, median)
    else median(x))
  out$StErr <- lapply(out$sims.list, function(x) {
    f <- function(x) {
      acf.x <- drop(acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
      acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
      ess <- length(x)/(1 + 2 * sum(acf.x))
      sqrt(var(x)/ess)
    }
    if (is.matrix(x))
      apply(x, 2, f)
    else f(x)
  })
  out$StDev <- lapply(out$sims.list, function(x) if (is.matrix(x))
    apply(x, 2, sd)
    else sd(x))
  # idem for mode
  out$modes <- lapply(out$sims.list, function(x) {
    m <- function(x) {
      d <- density(x, bw = "nrd", adjust = 3, n = 1000)
      d$x[which.max(d$y)]
    }
    if (is.matrix(x))
      apply(x, 2, m)
    else m(x)
  })


  if(!RE_ind && jags.data$ncU>1){
      Sigma2 <- matrix(0, ncol = jags.data$ncU, nrow = jags.data$ncU)
      Sigma2[lower.tri(Sigma2, TRUE)] <- out$coefficients$Sigma2
      Sigma2 <- Sigma2 + t(Sigma2)
      diag(Sigma2) <- diag(Sigma2)/2
      out$coefficients$Sigma2 <- Sigma2
      dimnames(out$coefficients$Sigma2) <- list(colnames(U),colnames(U))
      # idem for mode
      Sigma2 <- matrix(0, ncol = jags.data$ncU, nrow = jags.data$ncU)
      Sigma2[lower.tri(Sigma2, TRUE)] <- out$modes$Sigma2
      Sigma2 <- Sigma2 + t(Sigma2)
      diag(Sigma2) <- diag(Sigma2)/2
      out$modes$Sigma2 <- Sigma2
      dimnames(out$modes$Sigma2) <- list(colnames(U),colnames(U))
  }else{
    if(jags.data$ncU>1){
      out$coefficients$Sigma2 <- diag(out$coefficients$Sigma2)
      out$modes$Sigma2 <- diag(out$modes$Sigma2)
      dimnames(out$coefficients$Sigma2) <- dimnames(out$modes$Sigma2) <- list(colnames(U),colnames(U))
    }else{
      names(out$coefficients$Sigma2) <- names(out$modes$Sigma2) <- colnames(U)
    }
  }

  #---- End of the function defining the class and retruning the output
  class(out) <- "BQt"
  out

  }
