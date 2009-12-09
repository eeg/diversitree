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

# Rich's fix for NaN problem, merged in from diversitree.dev (11/19/09)
# problems with small D's on very short branches still remain and appear as a too-deep recursion error
#eps <- RTOL <- ATOL <- 1e-8
RTOL <- ATOL <- 1e-8
eps <- 0        # changed 12/7/09; this is how 0.4-1 has it
solve.gse2.C.guts <- function(y, len, pars, t0)
{
    ret <- t(gse2.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])

    if ( all(ret[,4:6] > eps) )
    {
        # should compare all of 4, 5, and 6 on next line
        q <- ret[cbind(seq_along(len), as.integer(ret[,5] > ret[,6]) + 5)]
        ret[,4:6] <- ret[,4:6] / q
        cbind(log(q), ret, deparse.level=0)
    } else
    {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)
        ret1 <- Recall(y, len1, pars, t0)
        ret2 <- Recall(ret1[n1,2:7], len2, pars, t0 + ti)
        ret2[,1] <- ret2[,1] + ret1[n1,1]
        rbind(ret1[-n1,], ret2)
    }
}

# a wrapper to reshape the answer from bisse.branches into old solve form
solve.gse2.C <- function(y, len, pars)
{
    ans <- solve.gse2.C.guts(y, len, pars, 0)
    ans[,5:7] <- ans[,5:7] * exp(ans[,1])
    ans <- t(ans)
    ans[1,] <- len
    ans
}

solve.gse2 <- solve.gse2.C
