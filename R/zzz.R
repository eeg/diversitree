.onLoad <- function(...) {
  assign("bisse.ode",
         make.ode("derivs", "diversitreeGSE", "initmod", 4),
         asNamespace("diversitreeGSE"))

  assign("gse2.ode",
         make.ode("gse2_derivs", "diversitreeGSE", "gse2_initmod", 6),
         asNamespace("diversitreeGSE"))
}
