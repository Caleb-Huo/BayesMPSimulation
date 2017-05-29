##' @export
myLimma <- function(dataExp, labels, seed = 15213) {
    library(limma)
    set.seed(seed)
    X <- matrix(, nrow = nrow(dataExp[[1]]), ncol = length(dataExp))
    for (i in 1:ncol(X)) {
        aData = dataExp[[i]]
        alabel = labels[[i]]
        labelData <- numeric(sum(sapply(alabel, length)))
        labelData[alabel[[1]]] = 0
        labelData[alabel[[2]]] = 1
        design = cbind(rep(1, length(labelData)), labelData)
        fit <- lmFit(aData, design)
        fit <- eBayes(fit)
        Amean <- fit$coefficients[, 2]
        X[, i] <- 0
        X[Amean > 0, i] <- -qnorm(fit$p.value[Amean > 0, 2]/2)
        X[Amean < 0, i] <- qnorm(fit$p.value[Amean < 0, 2]/2)
    }
    X
}
