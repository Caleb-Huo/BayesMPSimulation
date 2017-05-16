##' @export
simu.data.randomeffect_varyN2 <- function(K = 10, N1 = rep(50, K), N2 = rep(50, K), G = 10000, opdirectionRate = 0.01, 
    g = round(G * 0.1), p = rep(1/K, K), sigma = 1, sigmaRandom = 0.01, rho = 0.5, clust.size = 20, n.clust = 200, 
    rho.prior = 0.5, df.prior = clust.size * 3) {

	labels = NULL		
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
		
        controlLabel = 1:N1[i]
        caseLabel = (1:N2[i]) + N1[i]        
        labels[[i]] <- list(controlLabel = controlLabel, caseLabel = caseLabel)
		
    }
    truth = truthFixWithDirection
    truth[!truthSelection] = 0
    return(list(data = result, truth = truth, labels = labels))
}

