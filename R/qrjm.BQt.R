#' \code{qrjm.BQt} fits quantile regression joint model
#'
#' Function using JAGS to estimate the quantile regression joint model assuming asymmetric Laplace distribution for residual error.
#' Joint modeling concers longitudinal data and time-to-event
#'
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param formSurv survival formula as formula in survival package for latency submodel
#' @param survMod specifying the baseline risk function for Cox proportional hazard model (only "weibull" is available until now)
#' @param param shared association including in joint modeling: the classical shared random effects or the current value denoting by "sharedRE" (default) or "value", respectively.
#' @param timeVar string specify the names of time variable (time of repeated measurements)
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param RE_ind Boolean denoting if the random effects are assumed independent ; default is \code{FALSE}
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is 5000
#' @param quiet see rjags package
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
#' @author Antoine Barbieri
#'
#' @import rjags lcmm coda joineR
#'
#' @export
#'
#' @references Ming Yang, Sheng Luo, and Stacia DeSantis (2019).
#' \emph{Bayesian quantile regression joint models: Inference and dynamic predictions}.
#' Statistical Methods in Medical Research, 28(8):2524-2537. doi: 10.1177/0962280218784757.
#'
#' @examples
#'
#' \dontrun{
#' #---- use the data 'aids' from joineR package
#' data("aids", package = "joineR")
#'
#' #---- Fit quantile regression joint model for the first quartile
#' qrjm_25 <- qrjm.BQt(formFixed = CD4 ~ obstime,
#'                     formRandom = ~ obstime,
#'                     formGroup = ~ id,
#'                     formSurv = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                     survMod = "weibull",
#'                     param = "value",
#'                     timeVar= "obstime",
#'                     data = aids,
#'                     tau = 0.25)
#'
#' #---- Get the estimated coefficients
#' qrjm_25$coefficients
#'
#' #---- Summary of output
#' summary.BQt(qrjm_25)
#' }
#'
qrjm.BQt <- function(formFixed,
                     formRandom,
                     formGroup,
                     formSurv,
                     survMod = "weibull",
                     param = "value",
                     timeVar,
                     data,
                     tau,
                     RE_ind = FALSE,
                     n.chains = 1,
                     n.iter = 10000,
                     n.burnin = 5000,
                     n.thin = 5,
                     n.adapt = 5000,
                     quiet = FALSE,
                     C = 1000){

  # #
  # #   -- To do
  # #   verify with value.IG
  # #   propose shared random effects
  # #   add a stopping convergence criteria
  # #   initialize the values of parameter chaines
  # #
  #
  # data("aids", package = "joineR")
  #
  # formFixed = CD4 ~ obstime
  # formRandom = ~ obstime
  # formGroup = ~ id
  # formSurv = Surv(time, death) ~ drug + gender + prevOI + AZT
  # survMod = "weibull"
  # param = "sharedRE"
  # timeVar= "obstime"
  # data = aids
  # tau = 0.25
  # RE_ind = FALSE
  # n.chains = 1
  # n.iter = 2000
  # n.burnin = 1000
  # n.thin = 1
  # n.adapt = 5000
  # quiet = FALSE
  # C = 1000



  #-- data management

  # control
  lag = 0
  precision <- 1/100

  # longitudinal part
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  # use lcmm function to initiated values
  tmp_model <- lcmm::hlme(fixed = formFixed ,
                          random= formRandom,
                          subject = all.vars(formGroup),
                          data = data)
  # prior beta parameters
  priorMean.beta <- tmp_model$best[1:ncol(X)]
  priorTau.beta <- diag(rep(precision,length(priorMean.beta)))

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
                    priorA.sigma = precision,
                    priorB.sigma = precision
                    )

  if(jags.data$ncU==1)
    RE_ind <- TRUE
  if(RE_ind){
    jags.data <- c(jags.data,
                   list(priorA.Sigma2 = precision,
                        priorB.Sigma2 = precision
                   )
    )
  }else{
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(precision, ncol(U))),
                        priorK.Sigma2 = ncol(U),
                        mu0 = rep(0, ncol(U))
                   )
    )
  }

  # survival part
  tmp <- data[c(all.vars(formGroup),all.vars(formSurv))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(formSurv)][, 2]   # vector of event indicator (delta)
  nTime <- length(Time)                   # number of subject having Time
  zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- model.frame(formSurv, data = tmp)
  Z <- model.matrix(formSurv, mfZ)
  # fill the jags data
  priorMean.alpha <- rep(0, ncol(Z))
  priorTau.alpha <- diag(rep(precision, ncol(Z)))
  jags.data <- c(jags.data,
                 list(C = C,
                      zeros = numeric(nTime),
                      Time = Time,
                      event = event,
                      Z = Z,
                      ncZ = ncol(Z),
                      priorMean.alpha = priorMean.alpha,
                      priorTau.alpha = priorTau.alpha,
                      priorTau.alphaA = precision)
                 )

  #--- shared current value case
  data.id <- data_long[!duplicated(id), ]
  if (!timeVar %in% names(data_long))
    stop("\n'timeVar' does not correspond to one of the columns in formulas")
  if (param %in% c("value")) {
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    mfX.id <- model.frame(formFixed, data = data.id)
    Xtime <- model.matrix(formFixed, mfX.id)
    mfU.id <- model.frame(formRandom, data = data.id)
    Utime <- model.matrix(formRandom, mfU.id)
    # if (one.RE)
    #   Utime <- cbind(Utime, rep(0, nrow(Utime)))
    jags.data <- c(jags.data, list(Xtime = Xtime, Utime = Utime))

    #-- approxitmation of the intergral via the Gaussian quadrature (Gauss Kronrod rule)
    gaussKronrod <-
      function (k = 15) {
        sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
        wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                  0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                  0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                  0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
        wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
                 0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
        if (k == 7)
          list(sk = sk[1:7], wk = wk7)
        else
          list(sk = sk, wk = wk15)
      }

    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    K <- length(sk)
    P <- Time/2
    st <- outer(P, sk + 1)
    id.GK <- rep(seq_along(Time), each = K)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- c(t(st))
    mfX <- model.frame(formFixed, data = data.id2)
    mfU <- model.frame(formRandom, data = data.id2)
    Xs <- model.matrix(formFixed, mfX)
    Us <- model.matrix(formRandom, mfU)

    jags.data <- c(jags.data, list(K = K, P = P, st = st, wk = wk, Xs = Xs, Us = Us))
  }


  #---- Model for JAGS

  # model to consider
  model <- switch(paste(survMod, RE_ind, param, sep = "/"),
                  # wishart distribution when RE are considered dependent
                  `weibull/FALSE/value` = jags_qrjm.weib.value,
                  # inverse gamma distribution when RE are considered independent
                  `weibull/TRUE/value` = jags_qrjm.weib.value.IG,
                  # wishart distribution when RE are considered dependent
                  `weibull/FALSE/sharedRE` = jags_qrjm.weib.sharedRE,
                  # inverse gamma distribution when RE are considered independent
                  `weibull/TRUE/sharedRE` = jags_qrjm.weib.sharedRE.IG
                  )

  # parameters to save in the sampling step
  parms_to_save <- c("alpha", "alpha.assoc", "beta", "sigma", "b", "prec.Sigma2")

  # complement given survMod
  if(survMod == "weibull"){
    jags.data <- c(jags.data, list(priorA.shape = precision, priorB.shape = precision))
    parms_to_save <- c(parms_to_save, "shape")
  }

  # # initial values to improve chain convergence
  # initial.values <- list(gamma = priorMean.gamma,
  #                        alpha = priorMean.alpha,
  #                        shape = 3)

  #-- write jags model in txt from R function
  working.directory = getwd()
  write.model.jags(model = model,
                   name_model = "jags_qrjm",
                   intitled = file.path(working.directory, "JagsModel.txt"),
                   Data = jags.data,
                   param = param)

  #---- use JAGS sampler
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

  #---- output building

  #-- MCMClist management

  #- arguments
  out <- list(data = data)
  out$control <- list(fit = fit,
                      formFixed,
                      formRandom,
                      formGroup,
                      tau = tau,
                      call_function = "qrjm.BQt",
                      I = I,
                      C = C,
                      param = param,
                      survMod = survMod)

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
                    alpha = sims.list$alpha,
                    shape = sims.list$shape,
                    prec.Sigma2 = sims.list$prec.Sigma2)
  colnames(sims.list$beta) <- colnames(X)
  colnames(sims.list$alpha)[1:length(colnames(Z))] <- colnames(Z)

  #- random effect management
  ind.bs <- grep("b[", colnames(sims.list$b), fixed = TRUE)
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
