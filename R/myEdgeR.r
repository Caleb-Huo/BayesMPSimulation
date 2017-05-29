##' @export
myEdgeR <- function(dataExp, labels, seed = 15213) {
    library(edgeR)
    set.seed(seed)
    X <- matrix(, nrow = nrow(dataExp[[1]]), ncol = length(dataExp))
    for (i in 1:ncol(X)) {
        adata = dataExp[[i]]
        alabel = labels[[i]]
        labelData <- numeric(sum(sapply(alabel, length)))
        labelData[alabel[[1]]] = 0
        labelData[alabel[[2]]] = 1
		group<-factor(labelData)
		
		d=DGEList(counts=adata,group=group)
		d <- calcNormFactors(d) # normalization part
		design=model.matrix(~group)
		d=estimateGLMCommonDisp(d)
		d=estimateGLMTrendedDisp(d)
		d=estimateGLMTagwiseDisp(d)
		fit <- glmFit(d, design)
		lrt <- glmLRT(fit, coef=2)
		p_value=lrt$table[,4]
		acoef <- fit$coefficient[,2]
	
		Z1 <- rep(1,length(p_value))
		Z1[acoef>0] = -qnorm(p_value[acoef>0]/2)
		Z1[acoef<0] = qnorm(p_value[acoef<0]/2)
		
        X[, i] <- Z1
    }
    X
}
