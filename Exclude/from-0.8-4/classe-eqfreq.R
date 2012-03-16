library(diversitree)

#--------------------------------------------------
# Define parameters
#--------------------------------------------------

pars.bn <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05, 0.2, 0.15, 0.22, 0.08)
names(pars.bn) <- diversitree:::argnames.bisseness(NULL)
# "lambda0" "lambda1" "mu0"     "mu1"     "q01"     "q10"     
# "p0c"     "p0a"     "p1c"     "p1a"    

# translate from bisseness to classe
bn.to.cl <- function(pars)
{
    pars.cl <- pars * NA
    names(pars.cl) <- diversitree:::argnames.classe(NULL, 2)
    # "lambda111" "lambda112" "lambda122" "lambda211" "lambda212" "lambda222"
    # "mu1"       "mu2"       "q12"       "q21"      

    pars.cl["lambda111"] <- pars["lambda0"] * (1 - pars["p0c"])
    pars.cl["lambda112"] <- pars["lambda0"] * pars["p0c"] * pars["p0a"]
    pars.cl["lambda122"] <- pars["lambda0"] * pars["p0c"] * (1 - pars["p0a"])
    pars.cl["lambda211"] <- pars["lambda1"] * pars["p1c"] * (1 - pars["p1a"])
    pars.cl["lambda212"] <- pars["lambda1"] * pars["p1c"] * pars["p1a"]
    pars.cl["lambda222"] <- pars["lambda1"] * (1 - pars["p1c"])
    pars.cl[c("mu1", "mu2", "q12", "q21")] <- pars[c("mu0", "mu1", "q01", "q10")]

    # check
    # sum(pars.cl[c("lambda111", "lambda112", "lambda122")]) == pars["lambda0"]
    # sum(pars.cl[c("lambda222", "lambda212", "lambda211")]) == pars["lambda1"]
    # pars.cl[7:10] == pars[3:6]

    return(pars.cl)
}

# translate from classe to bisseness
cl.to.bn <- function(pars)
{
    pars.bn <- rep(NA, 10)
    pars.bn[1] <- sum(pars[1:3])
    pars.bn[2] <- sum(pars[4:6])
    pars.bn[3:6] <- pars[7:10]
    pars.bn[7] <- sum(pars[2:3]) / pars.bn[1]
    pars.bn[8] <- pars[2] / sum(pars[2:3])
    pars.bn[9] <- sum(pars[4:5]) / pars.bn[2]
    pars.bn[10] <- pars[5] / sum(pars[4:5])

    return(pars.bn)
}

#--------------------------------------------------
# 1. Mooch stationary.freq.bisseness()
#--------------------------------------------------

diversitree:::stationary.freq.bisseness(pars.bn)
# 0.3501673

pars.cl <- bn.to.cl(pars.bn)
diversitree:::stationary.freq.bisseness(cl.to.bn(pars.cl))
# 0.3501673

#--------------------------------------------------
# 2. Solve dx/dt for classe
#--------------------------------------------------

# only for k = 2
stationary.freq.classe <- function(pars)
{
    g <- (sum(pars[1:3]) - pars[7]) - (sum(pars[4:6]) - pars[8])
    eps <- sum(pars[1:8]) * 1e-14
    ss1 <- pars[9]  + 2*pars[3] + pars[2]  # shift from 1
    ss2 <- pars[10] + 2*pars[4] + pars[5]  # shift from 2

    if ( abs(g) < eps )
    {
        if (ss1 + ss2 == 0) 
            0.5
        else
            ss2/(ss1 + ss2)
    } else {
        roots <- diversitree:::quadratic.roots(g, ss2 + ss1 - g, -ss2)
        roots <- roots[roots >= 0 & roots <= 1]
        if ( length(roots) > 1 )
            NA
        else
            roots
    }
}

stationary.freq.classe(pars.cl)
# 0.3501673

#--------------------------------------------------
# 3. Eigenvector for classe
#--------------------------------------------------

stationary.freq.classe.k <- function(pars, k)
{
    # diversitree:::check.pars.classe(pars, k)
    nsum <- k*(k+1)/2
    kseq <- seq_len(k)
    pars.lam <- pars[seq(1, nsum*k)]
    pars.mu <- pars[seq(nsum*k+1, (nsum+1)*k)]
    pars.q <- pars[seq((nsum+1)*k+1, length(pars))]

    ## will be the transition matrix
    A <- matrix(0, nrow=k, ncol=k)

    ## array indices of lambda's in parameter vector
    idx.lam <- cbind(rep(kseq, each=nsum), rep(rep(kseq, times=seq(k,1,-1)), k),
                 unlist(lapply(kseq, function(i) i:k)))
    ## transpose of matrix indices of q's in parameter vector
    idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), rep(kseq, each=k-1))

    ## take care of off-diagonal elements
    for (n in seq_len(nsum*k))
    {
        ## add this lambda to A[daughter states, parent state]
        ## (separate steps in case the daughter states are the same)
        r <- idx.lam[n,]
        A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
        A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
    }
    A[idx.q] <- A[idx.q] + pars.q

    ## fix the diagonal elements
    diag(A) <- 0
    diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) 
                   sum(pars.lam[seq((i-1)*nsum+1, i*nsum)]) - pars.mu[i]))

    ## continuous time, so the dominant eigenvalue is the largest one
    ## return its eigenvector, normalized
    evA <- eigen(A)
    i <- which(evA$values == max(evA$values))
    evA$vectors[,i] / sum(evA$vectors[,i])
}

stationary.freq.classe.k(pars.cl, 2)

#--------------------------------------------------
# More tests
#--------------------------------------------------

# special cases

pars.bn <- c(1.4, 1.4, 0.6, 0.6, 0.8, 0.8, 0.4, 0.2, 0.4, 0.2)  # 0.5
names(pars.bn) <- diversitree:::argnames.bisseness(NULL)
diversitree:::stationary.freq.bisseness(pars.bn)
stationary.freq.classe(bn.to.cl(pars.bn))

pars.cl <- c(1.2, 0.1, 0.1, 0.1, 0.1, 1.2, 0.6, 0.6, 0.8, 0.8) # 0.5
pars.cl <- c(1.0, 0.2, 0.2, 0.1, 0.1, 1.2, 0.6, 0.6, 0.8, 0.8) # 0.44
diversitree:::stationary.freq.bisseness(cl.to.bn(pars.cl))
stationary.freq.classe(pars.cl)
stationary.freq.classe.k(pars.cl, 2)

# enough runs to compare elapsed times

pars.many <- matrix(runif(1000*10, min=0, max=3), ncol=10)

system.time( for (i in seq_len(nrow(pars.many))) stationary.freq.classe(pars.many[i,]) )
#  user  system elapsed 
# 0.060   0.000   0.061 

system.time( for (i in seq_len(nrow(pars.many))) stationary.freq.classe.k(pars.many[i,], 2) )
#  user  system elapsed 
# 0.850   0.000   0.845 

system.time( for (i in seq_len(nrow(pars.many))) diversitree:::stationary.freq.bisseness(cl.to.bn(pars.many[i,])) )
#  user  system elapsed 
# 0.110   0.000   0.106 

### For k = 3, compare with geosse

pars.ge <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars.ge) <- diversitree:::argnames.geosse(NULL)
pars.cl <- rep(0, 27)
names(pars.cl) <- diversitree:::argnames.classe(NULL, 3)

pars.cl['lambda222'] <- pars.cl['lambda112'] <- pars.ge['sA']
pars.cl['lambda333'] <- pars.cl['lambda113'] <- pars.ge['sB']
pars.cl['lambda123'] <-  pars.ge['sAB']
pars.cl['mu2'] <- pars.cl['q13'] <- pars.ge['xA']
pars.cl['mu3'] <- pars.cl['q12'] <- pars.ge['xB']
pars.cl['q21'] <- pars.ge['dA']
pars.cl['q31'] <- pars.ge['dB']

diversitree:::stationary.freq.geosse(pars.ge)
stationary.freq.classe.k(pars.cl, 3)
# 0.2764899 0.4966617 0.2268484
