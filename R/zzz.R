.onLoad <- function(...) {
  assign("bisse.ode",
         make.ode("derivs", "diversitreeGSE", "initmod", 4, FALSE),
         asNamespace("diversitreeGSE"))

  assign("gse2.ode",
         make.ode("gse2_derivs", "diversitreeGSE", "gse2_initmod", 6, FALSE),
         asNamespace("diversitreeGSE"))
}
