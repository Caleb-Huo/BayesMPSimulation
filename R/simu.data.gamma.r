##' @export
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

