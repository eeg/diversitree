## Constrain an argument to another argument
constrain.par <- function(f, rel) {
  f <- match.fun(f)
  rel.i <- resolve.constraint(rel)

#  if ( inherits(f, c("fixed", "constrained")) )
#    stop("Cannot (yet) constrain a constrained function")
#  else
    free <- is.na(rel)

  g <- function(x, ...)
    f(x[rel.i], ...)
  class(g) <- c(class(f), "constrained")
  g
}

## Fix an argument to a particular value.
fix.par <- function(f, rel) {
  f <- match.fun(f)

#  if ( inherits(f, c("fixed", "constrained")) )
#    ## free <- environment(f)$free & is.na(rel)
#    stop("Cannot (yet) constrain a constrained function")    
#  else
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


#-------------------------------
# EEG

# allow simultaneous constraints and fixing
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
    # are else's really needed with stop?

    free <- is.na(constrain.rel) & is.na(fix.rel)

    g <- function(x, ...)
        f(con.fix(x, constrain.rel, fix.rel), ...)
    class(g) <- c(class(f), "constrained")
    g
}

# doesn't work when constraining an unfixed parameter
con.fix <- function(x, constrain.rel, fix.rel)
{
    z = rep(NA, length(constrain.rel))

    i = which(!is.na(fix.rel))
    z[i] = fix.rel[i]

    i = which(!is.na(constrain.rel))
    z[i] = z[constrain.rel[i]]

    i = which(is.na(z))
    if (length(i) != length(x))
        stop("error in con.fix()")
    z[i] = x
    z
}

con.fix <- function(x, constrain.rel, fix.rel)
{
    z = rep(NA, length(constrain.rel))

    i = which(is.na(constrain.rel) & is.na(fix.rel))
    z[i] = x

    i = which(!is.na(fix.rel))
    z[i] = fix.rel[i]

    i = which(!is.na(constrain.rel))
    z[i] = z[constrain.rel[i]]

    return(z)
}



#--------------------------------------------------
# tests
#-------------------------------------------------- 

f = function(x) print(x)

x = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
f(x)

my.con = c(NA, 1, 1, NA, NA, 7, NA, NA)
f.con = constrain.par(f, my.con)
y = c(0.1, 0.4, 0.5, 0.7, 0.8)
f.con(y)      # 0.1 0.1 0.1 0.4 0.5 0.7 0.7 0.8

my.fix = c(0.01, NA, NA, NA, 0.05, NA, 0.07, NA)
f.fix = fix.par(f, my.fix)
y = c(0.2, 0.3, 0.4, 0.6, 0.8)
f.fix(y)      # 0.01 0.20 0.30 0.40 0.05 0.60 0.07 0.80

### constrain and fix
y = c(0.4, 0.8)

f.con.fix = constrain.par(f.fix, my.con)
f.con.fix(y)   # wrong: 0.01 0.40 0.40 0.40 0.05 0.80 0.07   NA

f.fix.con = fix.par(f.con, my.fix)
f.fix.con(y)   # wrong: 0.01 0.01 0.01 0.40 0.80   NA   NA 0.05


my.con = c(NA,    1,  1, NA,   NA,  7,   NA, NA)
my.fix = c(0.01, NA, NA, NA, 0.05, NA, 0.07, NA)
y = c(0.4, 0.8)
# desired:
f.con.fix = constrain.fix.par(f, my.con, my.fix)
f.con.fix(y)   # should be 0.01, 0.01, 0.01, 0.4, 0.05, 0.07, 0.07, 0.8

my.con = c(NA,    1,  1, NA,   NA,  7, NA, NA)
my.fix = c(0.01, NA, NA, NA, 0.05, NA, NA, NA)
y = c(0.4, 0.7, 0.8)
f.con.fix = constrain.fix.par(f, my.con, my.fix)
f.con.fix(y)   # should be 0.01, 0.01, 0.01, 0.4, 0.05, 0.7, 0.7, 0.8
