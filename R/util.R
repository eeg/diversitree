## Constrain an argument to another argument
constrain.par <- function(f, rel) {
  f <- match.fun(f)
  rel.i <- resolve.constraint(rel)

  if ( inherits(f, c("fixed", "constrained")) )
    stop("Cannot (yet) constrain a constrained function")
  else
    free <- is.na(rel)

  g <- function(x, ...)
    f(x[rel.i], ...)
  class(g) <- c(class(f), "constrained")
  g
}

## Fix an argument to a particular value.
fix.par <- function(f, rel) {
  f <- match.fun(f)

  if ( inherits(f, c("fixed", "constrained")) )
    ## free <- environment(f)$free & is.na(rel)
    stop("Cannot (yet) constrain a constrained function")    
  else
    free <- is.na(rel)

  g <- function(x, ...)
    f(spread(x, rel), ...)
  class(g) <- c(class(f), "fixed")
  g
}

spread <- function(x, rel) {
  idx <- is.na(rel)
  j <- which(!idx)
  y <- x[match(seq_along(rel), which(idx))]
  y[j] <- rel[j]
  y
}

resolve.constraint <- function(rel) {
  i <- is.na(rel)
  rel[i] <- seq_len(sum(i))
  rel[!i] <- rel[rel[!i]]
  rel
}


## Both constrain and fix arguments
constrain.fix.par <- function(f, constrain.rel, fix.rel)
{
    f <- match.fun(f)

    if (length(constrain.rel) != length(fix.rel))
        stop("Constrain and Fix vectors must be of equal length")
    if (any(!is.na(constrain.rel) & !is.na(fix.rel)))
        stop("Cannot both constrain and fix the same parameter")
    if (!all(constrain.rel %in% c(seq_along(constrain.rel), NA)))
        stop("Invalid constrain vector")
    if ( inherits(f, c("fixed", "constrained")) )
        stop("Cannot constrain an already-constrained function")
    # could also add a warning if constraining a parameter to itself
    # are else's really needed with stop?

    free <- is.na(constrain.rel) & is.na(fix.rel)

    g <- function(x, ...)
        f(con.fix(x, constrain.rel, fix.rel), ...)
    class(g) <- c(class(f), "constrained")
    g
}

# input:  x = vector of free parameter values
# output: z = vector of all parameter values, respecting constraints/fixing
con.fix <- function(x, constrain.rel, fix.rel)
{
    z = rep(NA, length(constrain.rel))

    i = which(is.na(constrain.rel) & is.na(fix.rel))
    z[i] = x   # works but sometimes gives a warning
#--------------------------------------------------
### first time through, x has length of all params
#     if (length(z) == length(x))
#     {
#         z[i] = x[i]
#     } else
#     {
#         z[i] = x
#     }
#-------------------------------------------------- 

    i = which(!is.na(fix.rel))
    z[i] = fix.rel[i]

    i = which(!is.na(constrain.rel))
    z[i] = z[constrain.rel[i]]

    return(z)
}


## Sensible starting point for optimisation
starting.point <- function(tree, q.div=5) {
 fit <- suppressWarnings(birthdeath(tree))
 r <- fit$para[2]
 e <- fit$para[1]
 p <- rep(c(r/(1-e), r*e/(1-e), r/q.div), each=2)
 names(p) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
 p
}

## Protect a function from failure.  This is a little like
## tryQuietly, but we will also ensure finiteness.
protect <- function(f, fail.value, finite=TRUE) {
  function(...) {
    ret <- trySilent(f(...))
    if ( inherits(ret, "try-error") ||
         (finite && (is.na(ret) || !is.finite(ret))) )
      fail.value
    else
      ret
  }
}

