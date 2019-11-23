replace.inprod <-
  function (body.model, name_model, Data) {
    mt <- deparse(body.model, width.cutoff = 200L)
    # regression from X
    ncX <- Data$ncX
    rplc <- paste(paste("beta[", 1:ncX, "] * X[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(beta[1:ncX], X[i, 1:ncX])", rplc, mt, fixed = TRUE)
    # regression on random effects U
    if(name_model != "jags_lqm"){
      ncU <- Data$ncU
      rplc <- paste(paste("b[i, ", 1:ncU, "] * U[j, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    # end
    c("model", mt)
  }
