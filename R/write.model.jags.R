write.model.jags <-
  function (model, name_model, intitled, Data) {
    model <- replace.inprod(body(model), name_model, Data)
    writeLines(model, intitled)
  }
