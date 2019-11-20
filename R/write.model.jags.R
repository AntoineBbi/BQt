write.model.jags <-
  function (model, intitled, Data) {
    model <- replace.inprod(body(model), Data)
    writeLines(model, intitled)
  }
