### make.ode() is in bisse-eqs.R

### needs changes
#--------------------------------------------------
# ## The R version of gse2-eqs.c:
# gse2.eqs <- function(t, x, pars)
# {
#     E.1  <- x[1]
#     E.2  <- x[2]
#     E.3  <- x[3]
#     D.N1 <- x[4]
#     D.N2 <- x[5]
#     D.N3 <- x[6]
# 
#     sA <- pars[1]
#     sB <- pars[2]
#     sAB
#     xA <- pars[3]
#     xB <- pars[4]
#     dA <- pars[5]
#     dB <- pars[6]
# 
#     need to include sAB
#     list(c(
#       -(sA+sB+xA+xB)*E.1 + xA* E.3 + xB* E.2 + sA*E.1*E.2 + sB*E.1*E.3,
#       -(sA + dA + xA) * E.2 + xA + dA * E.1 + sA * E.2 * E.2,
#       -(sB + dB + xB) * E.3 + xB + dB * E.1 + sB * E.3 * E.3,
#       -(sA+sB+xA+xB)*D.N1 + xA*D.N3 + xB*D.N2 + sA*(E.2*D.N1+E.1*D.N2) + sB*(E.3*D.N1 + E.1*D.N3),
#       -(sA + dA + xA) * D.N2 + dA * D.N1 + 2 * sA * D.N2 * E.2,
#       -(sB + dB + xB) * D.N3 + dB * D.N1 + 2 * sB * D.N3 * E.3
#     ))
# }
# 
# ## The Jacobian is not used
# #make.jacobian <- function(t, x, pars) { ... }
# 
# solve.R <- function(y, t, pars)
#   t(ode(y, c(0, t), gse2.eqs, pars))[-1,]
#-------------------------------------------------- 
# note that solve.C returns times as the first column and solve.R doesn't


### default tolerances specified in bisse-eqs.R
# RTOL <- 1e-8
# ATOL <- 1e-8

solve.gse2.C <- function(y, t, pars, rtol=RTOL, atol=ATOL)
  gse2.ode(y, c(0, t), pars, rtol=rtol, atol=atol)[,-1]

solve.gse2 <- solve.gse2.C
