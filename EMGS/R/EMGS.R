#' Perform EMGS method
#' 
#' This function performs EMGS method
#' 
#' @param X n by p matrix of data
#' @param v0 an increasing sequence of parameter v0, standard deviation of the spike distribution
#' @param v1 parameter v1, standard deviation of the slab distribution
#' @param lambda parameter for the diagonal elements
#' @param a prior for the sparsity parameter pi
#' @param b prior for the sparsity parameter pi
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param verbose logical indicator of printing information at each iteration
#' @param savepath logical indicator of saving the estimator at each iteration in the EM algorithm
#' @param copula logical indicator of copula model
#' @param copular_itr number of Monte Carlo iterations in simulating S
#' @param group grouping information
#' @param a_tau prior for the group scaling parameter
#' @param b_tau prior for the group scaling parameter
#' @param weights a vector of weights given to each observation
#' @return an fitted EMGS object
#' @references Li, Z. R., & McCormick, T. H. (2017). \emph{An Expectation Conditional Maximization approach for Gaussian graphical models}. arXiv preprint arXiv:1709.06970.
#' @examples
#' 
EMGS <- function(X, v0, v1, lambda, a, b, epsilon = 1e-5, maxitr = 1e4, verbose=TRUE, savepath = FALSE, copula = FALSE, copula_itr = 20, group=NULL, a_tau = 2, b_tau = 1, weights=NULL){
    
    if(copula){
        mranks <- rep(0, dim(X)[2])
        ranks <- X
        for(i in 1:dim(X)[2]){
            mranks[i] <- length(unique(X[, i])) 
            ranks[, i] <- match(X[, i], sort(unique(X[, i]))) - 1
        }
        latent_init <- qnorm(t(t(ranks+1) / (mranks+1)))
    }else{
        ranks <- X * 0
        mranks <- rep(0, dim(X)[2])
        latent_init = X
    }
    n <- dim(X)[1]
    if(is.null(weights)) weights <- rep(1/n, n)
    S <- t(X) %*% diag(weights) %*% X * n
    if(is.null(group)){
        group <- rep(0, dim(X)[2])
        exist.group <- 0
    }else{
        group <- match(group, unique(group)) - 1
        exist.group <- length(unique(group))
    }
    thin <- 1
    out <- .Call("_EMGS", as.matrix(X), S, v0, v1, lambda, a, b, epsilon, verbose,maxitr, savepath, copula, copula_itr, ranks, mranks, as.matrix(latent_init), thin, exist.group, group, a_tau, b_tau, PACKAGE = "EMGS")
    out$path <- out$EZ
	out$path[out$path >= 0.5] <- 1
	out$path[out$path < 0.5] <- 0
	for(i in 1:dim(out$EZ)[[3]]) diag(out$path[, , i]) <- 0
    out$v0 = v0
    out$v1 = v1
    out$lambda = lambda
    out$a = a
    out$b = b
    out$epsilon = epsilon
    out$maxitr = maxitr
    out$group = group

    if(exist.group > 0){
    taulist <- array(NA, dim = c(exist.group, exist.group, length(v0)))
    for(i in 1:length(v0)){
        for(j in 1:exist.group){
            for(k in 1:exist.group){
                taulist[j, k, i] <- out$tau[which(group == j-1)[1], which(group == k-1)[1], i]
            }
        }
    }
    out$tau_compact <- taulist
    }
    class(out) = "emgs"
    return(out)
}

#' Perform EMGS method with cross validation
#' 
#' This function performs EMGS method with cross validation
#' 
#' @param X n by p matrix of data
#' @param v0 an increasing sequence of parameter v0, standard deviation of the spike distribution
#' @param v1 parameter v1, standard deviation of the slab distribution
#' @param lambda parameter for the diagonal elements
#' @param a prior for the sparsity parameter pi
#' @param b prior for the sparsity parameter pi
#' @param K number of cross validation folds.
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param copula logical indicator of copula model
#' @param copular_itr number of Monte Carlo iterations in simulating S
#' @param group grouping information
#' @param ... other parameters passed to \code{\link{EMGS}}
#' @return  fitted EMGS objects with largest likelihood and corresponding to 1-SE rule.
#' @references Li, Z. R., & McCormick, T. H. (2017). \emph{An Expectation Conditional Maximization approach for Gaussian graphical models}. arXiv preprint arXiv:1709.06970.
#' @examples
#' 
cv.EMGS <- function(X, v0s, v1, lambda, a, b, epsilon = 1e-5, K=5, maxitr = 1e4,  copula, copula_itr, group=NULL, ...){
    getLL <- function(n, S, omega, copula){
        k <- dim(omega)[3]
        ll <- rep(NA, k)
        for(i in 1:k){
            ll[i] <- -n * log(det(omega[,,i])) + sum(diag(S%*%omega[,,i]))
        }            
        return(ll)
    }
    n <- dim(X)[1]
    X <- X[sample(n, n), ]
    ll <- matrix(0, length(v0s), K)
    for(k in 1:K){
        sub <- X[(1:n) %% K != (k-1), ]
        test <-  X[(1:n) %% K == (k-1), ]
        tmp <- EMGS(sub, v0s, v1, lambda, a, b, epsilon, maxitr, verbose=FALSE, savepath = FALSE, copula, copula_itr, group, ...)
        S <- t(test) %*% test
        ll[, k] <- getLL(n = dim(test)[1], S = S, omega = tmp$omega, copula) 
    }
    ll.mean <- apply(ll, 1, mean)
    ll.sd <- apply(ll, 1, sd) / sqrt(K - 1)
    v0.min <- v0s[which.min(ll.mean)]
    ll.1se <- ll.mean - (min(ll.mean) + ll.sd[which.min(ll.mean)])
    v0.1se<- v0s[max((1:length(v0s))[ll.1se <= 0])]

    fit.min <- EMGS(X=X, v0 = v0.min, v1=v1, lambda=lambda, a=a, b=b, epsilon=epsilon, maxitr=maxitr, verbose=FALSE, group = group, ...)
    fit.1se <- EMGS(X=X, v0 = v0.1se, v1=v1, lambda=lambda, a=a, b=b, epsilon=epsilon, maxitr=maxitr, verbose=FALSE, group = group, ...)
    return(list(ll = ll, v0s = v0s, v0.min = v0.min, v0.1se = v0.1se, fit.min = fit.min, fit.1se = fit.1se))
}


#' Perform Graphical lasso with cross validation
#' 
#' This function performs graphical lasso with cross validation
#' 
#' @param X n by p matrix of data
#' @param rholist list of tunning parameters
#' @param K number of cross validation folds.
#' @return  fitted objects with largest likelihood and corresponding to 1-SE rule.
#' @examples
#'
cv.glasso <- function(X, rholist = NULL, K = 5){
    require("glasso")
    getLL <- function(X, omega){
        n <- dim(X)[1]
        k <- length(omega) 
        ll <- rep(NA, k)
        for(i in 1:k){
            ll[i] <- -n * log(det(omega[[i]])) + sum(diag(t(X)%*%X%*%as.matrix(omega[[i]])))
        }
        return(ll)
    }
    X <- scale(X)
    S <- cor(X)
    if(is.null(rholist)) {
        lambda.min.ratio = 0.01
        d <- dim(X)[2]
        lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
        lambda.min = lambda.min.ratio * lambda.max
        rholist = rev(exp(seq(log(lambda.max), log(lambda.min), length = 30)))
    }
    n <- dim(X)[1]
    X <- X[sample(n, n), ]
    ll <- matrix(0, length(rholist), K)
    for(k in 1:K){
        sub <- X[(1:n) %% K != (k-1), ]
        test <-  X[(1:n) %% K == (k-1), ]
        tmp <- huge.glasso(sub, lambda = rholist)
        score <- getLL(X = test, omega = as.matrix(tmp$icov)) 
        ll[,k] <- score
    }
    ll.mean <- apply(ll, 1, mean)
    ll.sd <- apply(ll, 1, sd) / sqrt(K - 1)
    rho.min <- rholist[which.min(ll.mean)]
    ll.1se <- ll.mean - (min(ll.mean) + ll.sd[which.min(ll.mean)])
    rho.1se<- rholist[max((1:length(rholist))[ll.1se <= 0])]

    fit.min <-  huge.glasso(sub, lambda = rho.min) 
    fit.1se <- huge.glasso(sub, lambda = rho.1se)
    return(list(ll = ll, rholist = rholist, rho.min = rho.min, rho.1se = rho.1se, fit.min = fit.min, fit.1se = fit.1se))
}



  