#' \code{mqrjm.BQt} fits multiple quantile regression joint model
#'
#' Function using JAGS to estimate a multivariate quantile regression joint model assuming multivariate asymmetric Laplace
#' distribution for residual error. Specifically, the current value of median and the current value of range interquartile (Q3-Q1)
#' are considered as shared association between survival and longitudinal submodels.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param formSurv survival formula as formula in survival package for latency submodel
#' @param survMod specifying the baseline risk function for Cox proportional hazard model (only "weibull" is available until now)
#' @param param shared association including in joint modeling: only the current value denoting by "value" is considered.
#' When tau is a 3-length vector, the current range interquantile based on the two extrem values of tau and current value of the intermediate quantile are considered as shared association.
#' When tau is a 2-length vector, the current range interquantile is considered as shared association.
#' @param timeVar string specify the names of time variable (time of repeated measurements)
#' @param data dataset of observed variables
#' @param tau vector of quantiles to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped.
#' If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau},
#' otherwise the execution is stopped.
#' @param corr_structure string specifying the correlation structure between $Y|W_{tau}$: "free", "middle" or "none".
#' If 'none', the dependendance is only caught by the common covariance matrix of all random effects;
#' if 'middle', the model assumes a correlation structure based on both a common correlation parameter and the distances between considered quantiles's orders;
#' if 'free", a specific correlation parameter is assumed for each multivariate responses.
#' @param RE_ind if TRUE, the random effects are assumed to be independent between different quantile.
#' Otherwise (FALSE by default), correlation between quantile (markers) is also caught through common covariance matrix of random effects.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is NULL
#' @param precision variance by default for vague prior distribution
#' @param save_jagsUI If TRUE (by default), the output of jagsUI package is return by the function
#' @param parallel see jagsUI::jags() function
#' @param save_va If TRUE (is FALSE by default), the draws of auxilary variable W is returned by the function
#' @param C value used in the zero trick; default is 1000.
#'
#'
#' @return A \code{BQt} object is a list with the following elements:
#'  \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{random_effect}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific random effects}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#
#' @author Antoine Barbieri
#'
#' @import lqmm jagsUI survival
#'
#' @references Lea Petrella and Valentina Raponi (2019).
#' \emph{Joint estimation of conditional quantiles in multivariate linear regression models with an application to financial distress}.
#' Journal of Multivariate Analysis, 173:70-84. doi: 10.1016/j.jmva.2019.02.008.
#'
#' @references Elisabeth Waldmann and Thomas Kneib (2015).
#' \emph{Bayesian bivariate quantile regression}.
#' Statistical Modelling, 15(4):326-344. doi: 10.1177/1471082X14551247.
#'
#' @references Ming Yang, Sheng Luo, and Stacia DeSantis (2019).
#' \emph{Bayesian quantile regression joint models: Inference and dynamic predictions}.
#' Statistical Methods in Medical Research, 28(8):2524-2537. doi: 10.1177/0962280218784757.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' #---- Orthodont data from lqmm package
#' data("aids", package = "joineR")
#'
#' #---- Fit regression model for the first quartile
#' mqrjm_aids <- mqrjm.BQt(formFixed = CD4 ~ obstime,
#'                         formRandom = ~ obstime,
#'                         formGroup = ~ id,
#'                         formSurv = Surv(time, death) ~ drug + gender + prevOI + AZT,
#'                         survMod = "weibull",
#'                         param = "value",
#'                         timeVar= "obstime",
#'                         data = aids,
#'                         tau = c(0.25,0.5,0.75),
#'                         n.iter = 1000,
#'                         n.burnin = 500)
#'
#' #---- Get the estimated coefficients (posterior means)
#' mqrjm_aids$mean
#'
#' #---- Summary of output
#' summary.BQt(mqrjm_aids)
#'
#' #---- Traces of somes parameters
#' traceplot(x = mqrjm_aids$out_jagsUI, parameters = c("beta", "sigma", "alpha", "rho"))
#' }
mqrjm.BQt <- function(formFixed,
                      formRandom,
                      formGroup,
                      formSurv,
                      survMod = "weibull",
                      param = "value",
                      timeVar,
                      data,
                      tau,
                      corr_structure = "middle",
                      RE_ind = FALSE,
                      n.chains = 3,
                      n.iter = 10000,
                      n.burnin = 5000,
                      n.thin = 1,
                      n.adapt = NULL,
                      precision = 10,
                      save_jagsUI = TRUE,
                      parallel = FALSE,
                      save_va = FALSE,
                      C = 1000){

  # control
  lag = 0

  # Stopping condition
  if(prod(tau<1 && tau>0)!=1)
    stop("'tau' have to be a numerical vector with elements included between 0 and 1.\n")
  if(length(tau)>3)
    stop("'tau' have to be a numerical vector with length lower than 4.\n")
  if(sum(duplicated(tau))>0)
    stop("'Tau' have to be a numerical vector with no duplicated quantiles.\n")
  # if(RE_ind && corr_structure=="none"){
  #   cat("At least one correlation structure has to be defined, then 'RE_ind=FALSE'.\n")
  #   RE_ind <- FALSE
  # }
  if(corr_structure=="middle" && length(tau)==2){
    cat("The 'middle' correlation structure is equivalent to the 'free' one when two quantiles are considered (i.e. Q=2). Then, corr_structure=='free'.\n")
    corr_structure <- "free"
  }
  if(param=="sharedRE"){
    cat("Only param='value' is assumed for multiple quantile joint model.
        If two quantiles are considered, the latent strucutre is the interquantile range  ;
        if three quatiles are considered, the interquantile range between the two extreme quantiles and the current value of the middle one are considered as latent structure.\n")
    param <- "value"
  }

  ### data management

  #--- longitudinal part
  tau <- sort(tau)
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  ncU = ncol(U)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  Q <- length(tau)
  if(Q==2 && corr_structure=="middle")
    corr_structure = "free"
  # use lqmm function to initiated values
  cat("Initiation of parameter values using lqmm package. \n")
  tmp_model <- lqmm::lqmm(fixed = formFixed,
                          random = formRandom,
                          group = id,
                          tau = tau,
                          data = data_long)
  # prior beta parameters
  priorMean.beta <- t(coef(tmp_model))
  priorTau.beta <- diag(rep(1/10,length(priorMean.beta[1, ])))

  bis <- matrix(unlist(lqmm::ranef(tmp_model)), ncol = Q*ncU, byrow = F)
  bis[abs(bis)<.0001] <- 0
  init_sigma <- rep(1, Q)
  for(qq in 1:Q){
    init_sigma[qq] <- tmp_model[[qq]]$scale
  }
  initial.values <- list(beta = priorMean.beta,
                         sigma = init_sigma)
  if(RE_ind){
    initial.values$b1 <- bis[, 1:ncU]
    initial.values$b2 <- bis[, (ncU+1):(2*ncU)]
    if(Q==3)
      initial.values$b2 <- bis[, (ncU*2+1):(3*ncU)]
  }else{
    initial.values$b <- bis
  }

  # list of data jags
  jags.data <- list(y = matrix(rep(y,Q), ncol = Q, byrow = F),
                    X = X,
                    U = U,
                    tau = tau,
                    ncX = ncol(X),
                    ncU = ncol(U),
                    I = I,
                    Q = Q,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorA.sigma = 1/precision,
                    priorB.sigma = 1/precision
  )

  if(RE_ind){
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(1/precision, ncU)),
                        priorK.Sigma2 = ncU,
                        mu0 = rep(0, ncU)
                        )
                   )
  }else{
      jags.data <- c(jags.data,
                     list(priorR.Sigma2 = diag(rep(1/precision, ncU*Q)),
                          priorK.Sigma2 = ncU*Q,
                          mu0 = rep(0, ncU*Q)
                          )
                     )
  }

  prec.Sigma2 <- diag(1/unlist(lqmm::VarCorr(tmp_model)))
  if(RE_ind){
    initial.values$prec1.Sigma2 <- prec.Sigma2[1:ncU, 1:ncU]
    initial.values$prec1.Sigma2[initial.values$prec1.Sigma2 > 100] <- 100
    initial.values$prec2.Sigma2 <- prec.Sigma2[(ncU+1):(2*ncU), (ncU+1):(2*ncU)]
    initial.values$prec2.Sigma2[initial.values$prec2.Sigma2 > 100] <- 100
    if(Q==3){
      initial.values$prec3.Sigma2 <- prec.Sigma2[(ncU*2+1):(3*ncU), (ncU*2+1):(3*ncU)]
      initial.values$prec3.Sigma2[initial.values$prec3.Sigma2 > 100] <- 100
    }
  }else{
    initial.values$prec.Sigma2 <- prec.Sigma2
    initial.values$prec.Sigma2[initial.values$prec.Sigma2 > 100] <- 100
  }

  # manage the correlation between quantile
  if(corr_structure == "free"){
    if(Q==3){
      initial.values$rho12 <- 0.9
      initial.values$rho13 <- 0.9
      initial.values$rho23 <- 0.9
    }else{
      # Q=2
      initial.values$rho <- 0.9
    }
  }
  if(corr_structure == "middle"){
    jags.data$d12 <- tau[2]-tau[1]
    jags.data$d13 <- tau[3]-tau[1]
    jags.data$d23 <- tau[3]-tau[2]
  }


  #--- survival part
  tmp <- data[c(all.vars(formGroup),all.vars(formSurv))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(formSurv)][, 2]   # vector of event indicator (delta)
  nTime <- length(Time)                   # number of subject having Time
  zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- model.frame(formSurv, data = tmp)
  Z <- model.matrix(formSurv, mfZ)
  # use survival::coxph function to initiated values
  cat("> Initiation of survival parameter values using 'survival' package. \n")
  tmp_model <- survival::coxph(formSurv,
                               data = tmp,
                               x = TRUE)
  # Complete the jags data
  priorMean.alpha <- c(0, tmp_model$coefficients)
  priorTau.alpha <- diag(c(1/precision, 1/(precision*diag(tmp_model$var))))
  jags.data <- c(jags.data,
                 list(C = C,
                      zeros = numeric(nTime),
                      Time = Time,
                      event = event,
                      Z = Z,
                      ncZ = ncol(Z),
                      priorMean.alpha = priorMean.alpha,
                      priorTau.alpha = priorTau.alpha)
  )
  # initialisation values of survival parameters
  initial.values$alpha <- c(0, tmp_model$coefficients)
  # if(survMod=="weibull")
  #   initial.values$shape <- 1
  if(Q==2){
    initial.values$alpha.assoc <- 0
    jags.data <- c(jags.data,
                   list(priorMean.alpha.assoc = 0,
                   priorTau.alpha.assoc = 1/precision)
                   )
  }else{
    initial.values$alpha.assoc <- rep(0,Q-1)
    jags.data <- c(jags.data,
                   list(priorMean.alpha.assoc = rep(0,Q-1),
                   priorTau.alpha.assoc = diag(rep(1/precision,Q-1)))
                   )
  }

  #--- shared current value case
  data.id <- data_long[!duplicated(id), ]
  if (!timeVar %in% names(data_long))
    stop("'timeVar' does not correspond to one of the columns in formulas.\n")
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


  # write jags model (only Weibull baseline hazard function)
  model <- switch(paste(corr_structure, Q, RE_ind, sep = "/"),
                  `free/3/FALSE` = jags_3mqrjm_f.weib.value,
                  `free/3/TRUE`  = jags_3mqrjm_f_b.weib.value,
                  `free/2/FALSE` = jags_2mqrjm_f.weib.value,
                  `free/2/TRUE`  = jags_2mqrjm_f_b.weib.value,
                  `middle/3/FALSE` = jags_3mqrjm_m.weib.value,
                  `middle/3/TRUE`  = jags_3mqrjm_m_b.weib.value,
                  `none/3/FALSE` = jags_3mqrjm_n.weib.value,
                  `none/3/TRUE` = jags_3mqrjm_n_b.weib.value,
                  `none/2/FALSE` = jags_2mqrjm_n.weib.value
  )

  # parameters to save in the sampling step
  parms_to_save <- c("alpha", "alpha.assoc", "beta", "sigma")
  if(survMod == "weibull")
    parms_to_save <- c(parms_to_save, "shape")

  # specific parameters
  if(Q==3){
    if(corr_structure == "free"){
      parms_to_save <- c(parms_to_save, "rho12", "rho13", "rho23")
    }
    if(corr_structure == "middle"){
      jags.data$d12 <- tau[2]-tau[1]
      jags.data$d13 <- tau[3]-tau[1]
      jags.data$d23 <- tau[3]-tau[2]
      parms_to_save <- c(parms_to_save, "rho")
    }
    if(RE_ind){
      parms_to_save <- c(parms_to_save, "b1", "b2", "b3", "covariance.b1", "covariance.b2", "covariance.b3")
    }else{
      parms_to_save <- c(parms_to_save, "b", "covariance.b")
    }
  }else{# Q==2
      parms_to_save <- c(parms_to_save, "rho")
      if(RE_ind){
        parms_to_save <- c(parms_to_save, "b1", "b2", "covariance.b1", "covariance.b2")
      }else{
        parms_to_save <- c(parms_to_save, "b", "covariance.b")
      }
    }

  # useful to compute deviance in joint modeling
  if(save_va)
    parms_to_save <- c(parms_to_save, "W")

  # complement given survMod
  if(survMod == "weibull"){
    jags.data <- c(jags.data, list(priorA.shape = 1/precision, priorB.shape = 1/precision))
    parms_to_save <- c(parms_to_save, "shape")
  }

  #---- write jags model in txt from R function
  working.directory = getwd()
  write.model.jags(model = model,
                   name_model = "jags_mqrjm",
                   intitled = file.path(working.directory,"JagsModel.txt"),
                   Data = jags.data)

   # initalisation of chains
  if(n.chains==3)
    inits <- list(initial.values,
                  initial.values,
                  initial.values)
  if(n.chains==2)
    inits <- list(initial.values,
                  initial.values)
  if(n.chains==1)
    inits <- initial.values

  #---- JAGS sampler using jagsUI package
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = "JagsModel.txt",
                          inits = inits,
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          DIC = F)

  file.remove(file.path(working.directory, "JagsModel.txt"))

  #---- output building

  #-- MCMClist management

  #- arguments
  out <- list(data = data)
  out$control <- list(formFixed = formFixed,
                      formRandom = formRandom,
                      formGroup = formGroup,
                      formSurv = formSurv,
                      timeVar = timeVar,
                      tau = tau,
                      call_function = "mqrjm.BQt",
                      I = I,
                      C = C,
                      param = param,
                      survMod = survMod,
                      corr_structure = corr_structure,
                      n.chains = n.chains,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      parallel = parallel,
                      RE_ind = RE_ind,
                      event = event,
                      Time = Time)

  #- other outputs

  # sims.list output
  out$sims.list <- out_jags$sims.list
  if(!is.null(out$sims.list$out))
    out$sims.list$out <- NULL
  if(!RE_ind){
    out$sims.list$b <- NULL
  }else{
    out$sims.list$b1 <- out$sims.list$b2 <- out$sims.list$b3 <- NULL
  }

  # random effect output
  random_effect <- vector("list", Q)
  if(Q==3){
    names(random_effect) <- c(paste("tau", 100*tau[1], sep=""),
                              paste("tau", 100*tau[2], sep=""),
                              paste("tau", 100*tau[3], sep=""))
  }else{
    names(random_effect) <- c(paste("tau", 100*tau[1], sep=""),
                              paste("tau", 100*tau[2], sep=""))
  }

  if(!RE_ind){
    for(q in 1:Q){
      random_effect[[q]] <- list(postMeans = out_jags$mean$b[, ((q-1)*ncU+1):(q*ncU)],
                                 postSd = out_jags$sd$b[, ((q-1)*ncU+1):(q*ncU)])
      colnames(random_effect[[q]]$postMeans) <- colnames(U)
      colnames(random_effect[[q]]$postSd) <- colnames(U)
    }
  }else{
    random_effect[[1]] <- list(postMeans = out_jags$mean$b1, postSd = out_jags$sd$b1)
    random_effect[[2]] <- list(postMeans = out_jags$mean$b2, postSd = out_jags$sd$b2)
    colnames(random_effect[[1]]$postMeans) <-
      colnames(random_effect[[2]]$postMeans) <- colnames(U)
    colnames(random_effect[[1]]$postSd) <-
      colnames(random_effect[[2]]$postSd) <- colnames(U)
    if(Q==3){
      random_effect[[3]] <- list(postMeans = out_jags$mean$b3, postSd = out_jags$sd$b3)
      colnames(random_effect[[3]]$postMeans) <- colnames(random_effect[[3]]$postSd) <- colnames(U)
    }
  }
  out$random_effect <- random_effect

  # median : Posterior median of parameters (if mean, you can use mean instead of q50)
  out$median <- out_jags$q50
  if(!RE_ind){
    out$median$b <- NULL
  }else{
    out$median$b1 <- out$median$b2 <- NULL
    if(Q==3)
      out$median$b3 <- NULL
  }

  # mean : Posterior mean of parameters (if mean, you can use mean instead of q50)
  out$mean <- out_jags$mean
  if(!RE_ind){
    out$mean$b <- NULL
  }else{
    out$mean$b1 <- out$mean$b2 <- NULL
    if(Q==3)
      out$mean$b3 <- NULL
  }

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
  if(!RE_ind){
    out$StDev$b <- NULL
  }else{
    out$StDev$b1 <- out$StDev$b2 <- NULL
    if(Q==3)
      out$StDev$b3 <- NULL
  }

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  if(!RE_ind){
    out$Rhat$b <- NULL
  }else{
    out$Rhat$b1 <- out$Rhat$b2 <- NULL
    if(Q==3)
      out$Rhat$b3 <- NULL
  }
  # names
  names(out$mean$alpha) <-
    names(out$median$alpha) <-
    names(out$modes$alpha) <-
    names(out$StErr$alpha) <-
    names(out$Rhat$alpha) <-
    names(out$StDev$alpha) <- colnames(Z)

  colnames(out$mean$beta) <-
    colnames(out$median$beta) <-
    colnames(out$modes$beta) <-
    colnames(out$StErr$beta) <-
    colnames(out$Rhat$beta) <-
    colnames(out$StDev$beta) <- colnames(X)
  rownames(out$mean$beta) <-
    rownames(out$median$beta) <-
    rownames(out$modes$beta) <-
    rownames(out$StErr$beta) <-
    rownames(out$Rhat$beta) <-
    rownames(out$StDev$beta) <- paste("tau", as.character(tau*100), sep = "")
  rownames(out$mean$sigma) <-
    rownames(out$median$sigma) <-
    rownames(out$modes$sigma) <-
    rownames(out$StErr$sigma) <-
    rownames(out$Rhat$sigma) <-
    rownames(out$StDev$sigma) <- paste("tau", as.character(tau*100), sep = "")

  if(RE_ind){
    colnames(out$mean$covariance.b1) <-
      rownames(out$mean$covariance.b1) <-
      colnames(out$median$covariance.b1) <-
      rownames(out$median$covariance.b1) <-
      colnames(out$modes$covariance.b1) <-
      rownames(out$modes$covariance.b1) <-
      colnames(out$StErr$covariance.b1) <-
      rownames(out$StErr$covariance.b1) <-
      colnames(out$Rhat$covariance.b1) <-
      rownames(out$Rhat$covariance.b1) <-
      colnames(out$StDev$covariance.b1) <-
      rownames(out$StDev$covariance.b1) <-
      colnames(out$mean$covariance.b2) <-
      rownames(out$mean$covariance.b2) <-
      colnames(out$median$covariance.b2) <-
      rownames(out$median$covariance.b2) <-
      colnames(out$modes$covariance.b2) <-
      rownames(out$modes$covariance.b2) <-
      colnames(out$StErr$covariance.b2) <-
      rownames(out$StErr$covariance.b2) <-
      colnames(out$Rhat$covariance.b2) <-
      rownames(out$Rhat$covariance.b2) <-
      colnames(out$StDev$covariance.b2) <-
      rownames(out$StDev$covariance.b2) <-  colnames(U)
    if(Q==3){
      colnames(out$mean$covariance.b3) <-
        rownames(out$mean$covariance.b3) <-
        colnames(out$median$covariance.b3) <-
        rownames(out$median$covariance.b3) <-
        colnames(out$modes$covariance.b3) <-
        rownames(out$modes$covariance.b3) <-
        colnames(out$StErr$covariance.b3) <-
        rownames(out$StErr$covariance.b3) <-
        colnames(out$Rhat$covariance.b3) <-
        rownames(out$Rhat$covariance.b3) <-
        colnames(out$StDev$covariance.b3) <-
        rownames(out$StDev$covariance.b3) <- colnames(U)
    }
  }else{
    colnames(out$mean$covariance.b) <-
      rownames(out$mean$covariance.b) <-
      colnames(out$median$covariance.b) <-
      rownames(out$median$covariance.b) <-
      colnames(out$modes$covariance.b) <-
      rownames(out$modes$covariance.b) <-
      colnames(out$StErr$covariance.b) <-
      rownames(out$StErr$covariance.b) <-
      colnames(out$Rhat$covariance.b) <-
      rownames(out$Rhat$covariance.b) <-
      colnames(out$StDev$covariance.b) <-
      rownames(out$StDev$covariance.b) <- paste(rep(paste("tau", as.character(tau*100), sep = ""),
                                                    each = ncU),
                                                rep(colnames(U), Q),
                                                sep = ".")
  }

  # credible intervalles

  # survival part
  out$CIs$shape <- c(out_jags$q2.5$shape,
                     out_jags$q97.5$shape)
  names(out$CIs$shape) <- c("2.5%", "97.5%")

  out$CIs$alpha <- cbind(as.vector(t(out_jags$q2.5$alpha)),
                         as.vector(t(out_jags$q97.5$alpha)))
  rownames(out$CIs$alpha) <- colnames(Z)
  colnames(out$CIs$alpha) <- c("2.5%", "97.5%")

  if(Q==3){
    out$CIs$alpha.assoc <- cbind(as.vector(t(out_jags$q2.5$alpha.assoc)),
                                 as.vector(t(out_jags$q97.5$alpha.assoc)))
    rownames(out$CIs$alpha.assoc) <- c("current.Q", "current.rangeIQ")
    colnames(out$CIs$alpha.assoc) <- c("2.5%", "97.5%")
  }else{ # length(tau)==2
    out$CIs$alpha.assoc <- c(out_jags$q2.5$alpha.assoc,
                             out_jags$q97.5$alpha.assoc)
    names(out$CIs$alpha.assoc) <- c("2.5%", "97.5%")
  }

  # longitudinal part
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- paste(rep(paste("tau", as.character(tau*100), sep = ""),
                                      each = ncol(X)),
                                  rep(colnames(X), Q),
                                  sep = ".")
  out$CIs$sigma <- cbind(out_jags$q2.5$sigma,
                         out_jags$q97.5$sigma)
  rownames(out$CIs$sigma) <- paste("tau", as.character(tau*100), sep = "")
  if(corr_structure == "free"){
    if(Q==3){
      out$CIs$rho12 <- cbind(out_jags$q2.5$rho12,
                             out_jags$q97.5$rho12)
      out$CIs$rho13 <- cbind(out_jags$q2.5$rho13,
                             out_jags$q97.5$rho13)
      out$CIs$rho23 <- cbind(out_jags$q2.5$rho23,
                             out_jags$q97.5$rho23)
    }else{
      out$CIs$rho <- cbind(out_jags$q2.5$rho,
                           out_jags$q97.5$rho)
    }
  }
  if(corr_structure == "middle"){
    out$CIs$rho <- cbind(out_jags$q2.5$rho,
                         out_jags$q97.5$rho)
  }
  # only for diagonal elements of covariance matrix of random effects
  if(RE_ind){
    out$CIs$variances.b1 <- cbind(as.vector(diag(out_jags$q2.5$covariance.b1)),
                                  as.vector(diag(out_jags$q97.5$covariance.b1)))
    out$CIs$variances.b2 <- cbind(as.vector(diag(out_jags$q2.5$covariance.b2)),
                                  as.vector(diag(out_jags$q97.5$covariance.b2)))

    rownames(out$CIs$variances.b1) <- rownames(out$CIs$variances.b2) <- colnames(U)
    if(Q==3){
      out$CIs$variances.b3 <- cbind(as.vector(diag(out_jags$q2.5$covariance.b3)),
                                    as.vector(diag(out_jags$q97.5$covariance.b3)))
      rownames(out$CIs$variances.b3) <- colnames(U)
    }
  }else{
    out$CIs$variances.b <- cbind(as.vector(diag(out_jags$q2.5$covariance.b)),
                                 as.vector(diag(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$variances.b) <- paste(rep(paste("tau", as.character(tau*100), sep = ""),
                                               each = ncol(U)),
                                           rep(colnames(U), Q),
                                           sep = ".")
  }

  if(!RE_ind){
    if(corr_structure == "free"){
      if(Q==3){
        colnames(out$CIs$beta) <-
          colnames(out$CIs$sigma) <-
          colnames(out$CIs$rho12) <-
          colnames(out$CIs$rho13) <-
          colnames(out$CIs$rho23) <-
          colnames(out$CIs$variances.b) <- c("2.5%", "97.5%")
      }else{
        colnames(out$CIs$beta) <-
          colnames(out$CIs$sigma) <-
          colnames(out$CIs$rho) <-
          colnames(out$CIs$variances.b) <- c("2.5%", "97.5%")
      }
    }
    if(corr_structure == "middle"){
      colnames(out$CIs$beta) <-
        colnames(out$CIs$sigma) <-
        colnames(out$CIs$rho) <-
        colnames(out$CIs$variances.b) <- c("2.5%", "97.5%")
    }
    if(corr_structure == "none"){
      colnames(out$CIs$beta) <-
        colnames(out$CIs$sigma) <-
        colnames(out$CIs$variances.b) <- c("2.5%", "97.5%")
    }
  }else{
    if(corr_structure == "free"){
      if(Q==3){
        colnames(out$CIs$beta) <-
          colnames(out$CIs$sigma) <-
          colnames(out$CIs$rho12) <-
          colnames(out$CIs$rho13) <-
          colnames(out$CIs$rho23) <-
          colnames(out$CIs$variances.b1) <-
          colnames(out$CIs$variances.b2) <-
          colnames(out$CIs$variances.b3) <- c("2.5%", "97.5%")
      }else{# not implemented yet...
        colnames(out$CIs$beta) <-
          colnames(out$CIs$sigma) <-
          colnames(out$CIs$rho) <-
          colnames(out$CIs$variances.b1) <-
          colnames(out$CIs$variances.b2) <-
          colnames(out$CIs$variances.b3) <-  c("2.5%", "97.5%")
      }
    }
    if(corr_structure == "middle"){
      colnames(out$CIs$beta) <-
        colnames(out$CIs$sigma) <-
        colnames(out$CIs$rho) <-
        colnames(out$CIs$variances.b1) <-
        colnames(out$CIs$variances.b2) <-
        colnames(out$CIs$variances.b3) <- c("2.5%", "97.5%")
    }
    if(corr_structure == "none"){
      colnames(out$CIs$beta) <-
        colnames(out$CIs$sigma) <-
        colnames(out$CIs$variances.b1) <-
        colnames(out$CIs$variances.b2) <-
        colnames(out$CIs$variances.b3) <- c("2.5%", "97.5%")
    }
  }

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  #---- End of the function defining the class and retruning the output

  class(out) <- "BQt"
  out

}
