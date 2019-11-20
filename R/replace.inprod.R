replace.inprod <- 
  function (body.model, Data) {
    mt <- deparse(body.model, width.cutoff = 200L)
    # regression from X
    ncX <- Data$ncX
    rplc <- paste(paste("beta[", 1:ncX, "] * X[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(beta[1:ncX], X[i, 1:ncX])", rplc, mt, fixed = TRUE)
    c("model", mt)
  }
