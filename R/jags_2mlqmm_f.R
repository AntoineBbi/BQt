jags_2mlqmm_f <-
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
        V11[j] <- W[j, 1]*sigma[1]*c2[1]
        V22[j] <- W[j, 2]*sigma[2]*c2[2]
        V12[j] <- sqrt(V11[j]*V22[j])*rho
        # first quantile distribution
        y[j, 1] ~ dnorm(mu1[j], prec1[j])
        mu1[j] <- inprod(beta[1, 1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1[1]*W[j, 1]
        prec1[j] <- 1/V11[j]
        # conditional distrubtion for second quantile given y[j, 1]
        y[j, 2] ~ dnorm(mu2.knowing.y1[j], prec2[j])
        mu2[j] <- inprod(beta[2, 1:ncX], X[j, 1:ncX]) + inprod(b[i, (ncU+1):(2*ncU)], U[j, 1:ncU]) + c1[2]*W[j, 2]
        mu2.knowing.y1[j] <- mu2[j] + rho*sqrt(V22[j]/V11[j])*(y[j, 1]-mu1[j])
        prec2[j] <- 1/(V22[j]*(1-pow(rho, 2)) )
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
    rho ~ dunif(-1,1)
  }
