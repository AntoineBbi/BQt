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
        for(qq in 1:Q){
          # define object
          W[j, qq] ~ dexp(1/sigma[qq])
          prec[j, qq] <- 1/(W[j, qq]*sigma[qq]*c2[qq])
          # first quantile distribution
          y[j, 1] ~ dnorm(mu[j, qq], prec[j, qq])
          mu3[j] <- inprod(beta[3, 1:ncX], X[j, 1:ncX]) + inprod(b[i, (ncU*(qq-1)+1):(qq*ncU)], U[j, 1:ncU]) + c1[qq]*W[j, qq]
        }#end of qq loop
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
