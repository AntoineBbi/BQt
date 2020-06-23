jags_3mlqmm_n <-
  function (...) {
    for(q in 1:Q){
      c1[q] <- (1-2*tau[q])/(tau[q]*(1-tau[q]))
      c2[q] <- 2/(tau[q]*(1-tau[q]))
    }
    # likelihood
    for (i in 1:I){
      # longitudinal part
      for(j in offset[i]:(offset[i+1]-1)){
        # define object
        W[j, 1] ~ dexp(1/sigma[1])
        W[j, 2] ~ dexp(1/sigma[2])
        W[j, 3] ~ dexp(1/sigma[3])
        prec1[j] <- 1/(W[j, 1]*sigma[1]*c2[1])
        prec2[j] <- 1/(W[j, 2]*sigma[2]*c2[2])
        prec3[j] <- 1/(W[j, 3]*sigma[3]*c2[3])
        # first quantile distribution
        y[j, 1] ~ dnorm(mu1[j], prec1[j])
        mu1[j] <- inprod(beta[1, 1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1[1]*W[j, 1]
        # conditional distrubtion for second quantile given y[j, 1]
        y[j, 2] ~ dnorm(mu2[j], prec2[j])
        mu2[j] <- inprod(beta[2, 1:ncX], X[j, 1:ncX]) + inprod(b[i, (ncU+1):(2*ncU)], U[j, 1:ncU]) + c1[2]*W[j, 2]
        # # Conditional normal distribution for third quantile, conditional on both first and second quantile
        y[j, 3] ~ dnorm(mu3[j], prec3[j])
        mu3[j] <- inprod(beta[3, 1:ncX], X[j, 1:ncX]) + inprod(b[i, (ncU*2+1):(3*ncU)], U[j, 1:ncU]) + c1[3]*W[j, 3]
      }#end of j loop
      # random effects
      b[i, 1:(ncU*Q)] ~ dmnorm(mu0[], prec.Sigma2[, ])
    }#end of i loop
    # priors for parameters
    prec.Sigma2[1:(ncU*Q), 1:(ncU*Q)] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
    covariance.b <- inverse(prec.Sigma2[, ])
    for(qqq in 1:Q){
      beta[qqq, 1:ncX] ~ dmnorm(priorMean.beta[qqq, ], priorTau.beta[, ])
      sigma[qqq] ~ dgamma(priorA.sigma, priorB.sigma)
    }
  }
