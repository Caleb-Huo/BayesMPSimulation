################################### Functions # revision by Xiaoguang 2014/10/20, GSPH 309D guarantee the same effect size direction

library(MCMCpack)
library(truncnorm)


sample.correlated.data <- function(p, n, rho, df, sigma) {
    s <- matrix(rho, p, p)
    s <- s + diag(rep(1 - rho, p))
    cov.mat <- riwish(df, s)
    cor.mat <- cov2cor(cov.mat)
    eig <- eigen(cor.mat)
    L <- eig[[2]] %*% diag(sqrt(eig[[1]])) %*% t(eig[[2]])
    x <- matrix(rnorm(p * n, 0, sigma), p, n)
    xx <- L %*% x
    return(xx)
}

simu.data.correlated <- function(K = 10, N = 50, G = 10000, opdirectionRate = 0.01, g = round(G * 0.1), p = rep(1/K, 
    K), mu = matrix(runif(G * K, 0.5, 1) * sample(c(-1, 1), G, replace = TRUE), G, K) * sample(c(-1, 1), G * 
    K, replace = TRUE, prob = c(opdirectionRate, 1 - opdirectionRate)), sigma = 1, rho = 0.5, clust.size = 20, 
    n.clust = 200, rho.prior = 0.5, df.prior = clust.size * 3) {
    k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
    truth <- matrix(FALSE, G, K)
    for (i in 1:G) truth[i, sample(1:K, k[i])] <- TRUE
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
        mu.i[!truth[, i]] <- 0
        data0[, (1:N) + N] <- data0[, (1:N) + N] + mu.i
        result[[i]] <- data0
    }
    return(list(data = result, truth = truth))
}

simu.data.randomeffect <- function(K = 10, N = 50, G = 10000, opdirectionRate = 0.01, g = round(G * 0.1), p = rep(1/K, 
    K), sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, n.clust = 200, rho.prior = 0.5, df.prior = clust.size * 
    3) {
    truthFix = sample(c(-1, 1), G, replace = TRUE)
    # mufix = runif(G , a, b) * truthFix
    mufix = rtruncnorm(G, a = 0.5, b = Inf, mean = 1, sd = 1)
    muRandom = matrix(rnorm(G * K, mufix, sigmaRandom), ncol = K, nrow = G)
    truthDirection = sample(c(-1, 1), G * K, replace = TRUE, prob = c(opdirectionRate, 1 - opdirectionRate))
    truthFixWithDirection = matrix(truthFix * truthDirection, G, K)
    mu = muRandom * truthFixWithDirection
    
    k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
    truthSelection <- matrix(FALSE, G, K)
    for (i in 1:G) truthSelection[i, sample(1:K, k[i])] <- TRUE
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
        mu.i[!truthSelection[, i]] <- 0
        data0[, (1:N) + N] <- data0[, (1:N) + N] + mu.i
        result[[i]] <- data0
    }
    truth = truthFixWithDirection
    truth[!truthSelection] = 0
    return(list(data = result, truth = truth))
}


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

# generateS = simu.data.partial(K = 4, N = 20, G = 10000,typeselection=c(1),comb =
# c(0.04,0.02,0.02,0.2),sigma=1, sigmaRandom = 0.01) generateS = simu.data.partial(K = 4, N = 20, G =
# 10000,typeselection=c(1,2),firstK=2,comb = c(0.04,0.04,0.02,0.2),sigma=1, sigmaRandom = 0.01)

simu.data.gamma <- function(K = 10, N = 50, G = 10000, opdirectionRate = 0.01, g = round(G * 0.1), p = rep(1/K, 
    K), scale = 2, shape = 2, truncationGamma = 1, sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, 
    n.clust = 200, rho.prior = 0.5, df.prior = clust.size * 3) {
    truthFix = sample(c(-1, 1), G, replace = TRUE)
    # mufix = runif(G , a, b) * truthFix
    NN <- Gammad(scale = scale, shape = shape)
    NT <- Truncate(NN, lower = truncationGamma, upper = Inf)
    
    mufix = r(NT)(G)
    # rtruncnorm(G, a = 0.5, b = Inf, mean = 1, sd = 1)
    muRandom = matrix(rnorm(G * K, mufix, sigmaRandom), ncol = K, nrow = G)
    truthDirection = sample(c(-1, 1), G * K, replace = TRUE, prob = c(opdirectionRate, 1 - opdirectionRate))
    truthFixWithDirection = matrix(truthFix * truthDirection, G, K)
    mu = muRandom * truthFixWithDirection
    
    k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
    truthSelection <- matrix(FALSE, G, K)
    for (i in 1:G) truthSelection[i, sample(1:K, k[i])] <- TRUE
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
        mu.i[!truthSelection[, i]] <- 0
        data0[, (1:N) + N] <- data0[, (1:N) + N] + mu.i
        result[[i]] <- data0
    }
    truth = truthFixWithDirection
    truth[!truthSelection] = 0
    return(list(data = result, truth = truth))
}

## varying N for different studies
simu.data.randomeffect_varyN <- function(K = 10, N = rep(50, K), G = 10000, opdirectionRate = 0.01, g = round(G * 
    0.1), p = rep(1/K, K), sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, n.clust = 200, rho.prior = 0.5, 
    df.prior = clust.size * 3) {
    truthFix = sample(c(-1, 1), G, replace = TRUE)
    # mufix = runif(G , a, b) * truthFix
    mufix = rtruncnorm(G, a = 0.5, b = Inf, mean = 1, sd = 1)
    muRandom = matrix(rnorm(G * K, mufix, sigmaRandom), ncol = K, nrow = G)
    truthDirection = sample(c(-1, 1), G * K, replace = TRUE, prob = c(opdirectionRate, 1 - opdirectionRate))
    truthFixWithDirection = matrix(truthFix * truthDirection, G, K)
    mu = muRandom * truthFixWithDirection
    
    k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
    truthSelection <- matrix(FALSE, G, K)
    for (i in 1:G) truthSelection[i, sample(1:K, k[i])] <- TRUE
    ttt <- rep(0, G)
    ttt[1:(n.clust * clust.size)] <- rep(1:n.clust, clust.size)
    clust <- sample(ttt)
    result <- list()
    for (i in 1:K) {
        print(i)
        data0 <- matrix(rnorm(2 * G * N[i], 0, sigma), G, N[i] * 2)
        for (j in 1:n.clust) {
            data0[clust == j, ] <- sample.correlated.data(clust.size, N[i] * 2, rho.prior, df.prior, sigma)
        }
        mu.i <- mu[, i]
        mu.i[!truthSelection[, i]] <- 0
        data0[, (1:N[i]) + N[i]] <- data0[, (1:N[i]) + N[i]] + mu.i
        result[[i]] <- data0
    }
    truth = truthFixWithDirection
    truth[!truthSelection] = 0
    return(list(data = result, truth = truth))
}

## varying N within study.
simu.data.randomeffect_varyN2 <- function(K = 10, N1 = rep(50, K), N2 = rep(50, K), G = 10000, opdirectionRate = 0.01, 
    g = round(G * 0.1), p = rep(1/K, K), sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, n.clust = 200, 
    rho.prior = 0.5, df.prior = clust.size * 3) {
    truthFix = sample(c(-1, 1), G, replace = TRUE)
    # mufix = runif(G , a, b) * truthFix
    mufix = rtruncnorm(G, a = 0.5, b = Inf, mean = 1, sd = 1)
    muRandom = matrix(rnorm(G * K, mufix, sigmaRandom), ncol = K, nrow = G)
    truthDirection = sample(c(-1, 1), G * K, replace = TRUE, prob = c(opdirectionRate, 1 - opdirectionRate))
    truthFixWithDirection = matrix(truthFix * truthDirection, G, K)
    mu = muRandom * truthFixWithDirection
    
    k <- c(sort(sample(1:K, g, replace = TRUE), decreasing = TRUE), rep(0, G - g))
    truthSelection <- matrix(FALSE, G, K)
    for (i in 1:G) truthSelection[i, sample(1:K, k[i])] <- TRUE
    ttt <- rep(0, G)
    ttt[1:(n.clust * clust.size)] <- rep(1:n.clust, clust.size)
    clust <- sample(ttt)
    result <- list()
    for (i in 1:K) {
        print(i)
        data0 <- matrix(rnorm(G * (N1[i] + N2[i]), 0, sigma), G, (N1[i] + N2[i]))
        for (j in 1:n.clust) {
            data0[clust == j, ] <- sample.correlated.data(clust.size, (N1[i] + N2[i]), rho.prior, df.prior, 
                sigma)
        }
        mu.i <- mu[, i]
        mu.i[!truthSelection[, i]] <- 0
        data0[, (1:N2[i]) + N1[i]] <- data0[, (1:N2[i]) + N1[i]] + mu.i
        result[[i]] <- data0
    }
    truth = truthFixWithDirection
    truth[!truthSelection] = 0
    return(list(data = result, truth = truth))
}

