##' @export
simu.data.partial <- function(K = 10, N = 50, G = 10000, comb = c(0.02, 0.02, 0.02, 0.02), g = round(G * 0.1), 
    p = rep(1/K, K), typeselection = c(1, 2, 3, 4), firstK = K, sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, 
    n.clust = 200, rho.prior = 0.5, df.prior = clust.size * 3) {
    truthFix = sample(c(-1, 1), G, replace = TRUE)
    # mufix = runif(G , a, b) * truthFix
    mufix = rtruncnorm(G, a = 0.5, b = Inf, mean = 1, sd = 1)
    muRandom = matrix(rnorm(G * K, mufix, sigmaRandom), ncol = K, nrow = G)
    
    homoConcordant <- function(K, arow) {
        matrix(sample(c(-1, 1), arow, replace = TRUE) * rep(1, K * arow), nrow = arow, ncol = K)
    }
    # homoConcordant(5,5)
    
    studySpecific <- function(K, arow, afirstK = firstK) {
        ## firstK could be significant, say K=10 and firstK=4
        res = matrix(0, nrow = arow, ncol = K)
        for (i in 1:arow) {
            res[i, sample(1:afirstK, 1)] = sample(c(-1, 1), 1)
        }
        res
    }
    # studySpecific(5,5,3)
    
    partialStudySpecific <- function(K, arow) {
        res = matrix(0, nrow = arow, ncol = K)
        for (i in 1:arow) {
            res[i, 1:sample(1:(K - 1), 1)] = sample(c(-1, 1), 1)
        }
        res
    }
    # partialStudySpecific(5,5)
    
    discordant <- function(K, arow, rate = 0.4) {
        matrix(sample(c(1, 0), K * arow, replace = TRUE, prob = c(rate, 1 - rate)) * sample(c(1, -1), K * arow, 
            replace = TRUE), nrow = arow, ncol = K)
    }
    # discordant(5,5)
    
    functions <- c(homoConcordant, studySpecific, partialStudySpecific, discordant)[typeselection]
    
    DEtype = character(G)
    DEfactor = c("homoConcordant", "studySpecific", "partialStudySpecific", "discordant")[typeselection]
    pointer = 1
    
    truthFixWithDirection = matrix(0, G, K)
    
    for (i in 1:length(functions)) {
        anum = floor(comb[i] * G)
        anindex = pointer:(pointer + anum - 1)
        DEtype[anindex] = rep(DEfactor[i], anum)
        truthFixWithDirection[anindex, ] = functions[[i]](K, anum)
        pointer = pointer + anum
    }
    
    mu = muRandom * truthFixWithDirection
    
    ttt <- rep(0, G)
    ttt[1:(n.clust * clust.size)] <- rep(1:n.clust, clust.size)
    clust <- sample(ttt)
    result <- list()
    
    for (i in 1:K) {
        print(i)
        data0 <- matrix(rnorm(2 * G * N, 0, sigma), G, N * 2)
        for (j in 1:n.clust) {
            data0[clust == j, ] <- sample.correlated.data(clust.size, N * 2, rho.prior, df.prior, sigma)
        }
        mu.i <- mu[, i]
        data0[, (1:N) + N] <- data0[, (1:N) + N] + mu.i
        result[[i]] <- data0
    }
    truth = truthFixWithDirection
    
    return(list(data = result, truth = truth, DEtype = DEtype))
}

