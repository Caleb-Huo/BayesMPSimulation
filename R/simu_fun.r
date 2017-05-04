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

