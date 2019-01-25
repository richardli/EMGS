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
#' @param weights a vector of weights given to each observation, not used.
#' @param warm logical indicator for warm start
#' @return an fitted EMGS object
#' @references Li, Z. R., & McCormick, T. H. (2017). \emph{An Expectation Conditional Maximization approach for Gaussian graphical models}. arXiv preprint arXiv:1709.06970.
#' @examples
#' 
EMGS <- function(X, v0, v1, lambda, a, b, epsilon = 1e-5, maxitr = 1e4, verbose=TRUE, savepath = FALSE, copula = FALSE, copula_itr = 20, group=NULL, a_tau = 2, b_tau = 1, weights=NULL, warm = FALSE, Xtest = NULL){
    
    if(copula){
        mranks <- rep(0, dim(X)[2])
        ranks <- X
        for(i in 1:dim(X)[2]){
            mranks[i] <- length(unique(X[, i])) 
            ranks[, i] <- match(X[, i], sort(unique(X[, i]))) - 1
        }
        # latent_init <- qnorm(t(t(ranks+1) / (mranks+1)))
        latent_init <- huge::huge.npn(X)
        if(!is.null(Xtest)){
            mranks.test <- rep(0, dim(Xtest)[2])
            ranks.test <- Xtest
            for(i in 1:dim(Xtest)[2]){
                mranks.test[i] <- mranks[i]
                cuts <- sort(unique(X[, i]))
                for(ii in 1:dim(Xtest)[1]){
                    #    r1   r2    r3    r4
                    # ( rt1  )[ rt2)[ rt3)[rt4  )
                    ranks.test[ii, i] <- max(which(order(c(Xtest[ii, i], cuts))==1)-1, 0)
                }
            }
            latent_init.test <- huge::huge.npn(Xtest)
            Stest <- t(latent_init.test) %*% latent_init.test  
        }else{
            Xtest <- Stest <- ranks.test <- latent_init.test <- matrix(0, 1, 1)
            mranks.test <- rep(0, 1)
        }
    }else{
        ranks <- X * 0
        mranks <- rep(0, dim(X)[2])
        latent_init <- X
        ranks.test <- matrix(0, 1, 1)
        if(is.null(Xtest)){
            Xtest <- Stest <- latent_init.test <- matrix(0, 1, 1)
        }else{
            latent_init.test <- Xtest
            latent_init.test[is.na(latent_init.test)] <- 0
             Stest <- t(latent_init.test) %*% latent_init.test             
        }
        mranks.test <- rep(0, 1)
    }
    allmissing <- apply(X, 1, function(x) sum(is.na(x)))
    X <- X[which(allmissing < dim(X)[2]), ]
    n <- dim(X)[1]
    has.missing <- sum(is.na(X)) > 0
    if(is.null(weights)) weights <- rep(1/n, n)
    if(has.missing){
        tmp <- latent_init
        tmp[is.na(tmp)] <- 0
        S <- t(tmp) %*% diag(weights) %*% tmp * n
    }else{
        S <- t(latent_init) %*% diag(weights) %*% latent_init * n
    }
    if(is.null(group)){
        group <- rep(0, dim(X)[2])
        exist.group <- 0
    }else{
        group <- match(group, unique(group)) - 1
        exist.group <- length(unique(group))
    }
    thin <- 1
    out <- .Call("_EMGS", as.matrix(X), S, v0, v1, lambda, a, b, epsilon, verbose,maxitr, savepath, copula, copula_itr, ranks, mranks, as.matrix(latent_init), thin, exist.group, group, a_tau, b_tau, warm, Xtest, Stest, ranks.test, mranks.test, latent_init.test, has.missing, PACKAGE = "EMGS")
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
    if(!has.missing) out$X = NULL
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
    getLL <- function(n, cov, omega, copula){
        k <- dim(omega)[3]
        ll <- rep(NA, k)
        for(i in 1:k){
            # if(!copula){
                ll[i] <- - log(det(omega[,,i])) + sum(diag(cov%*%omega[,,i]))
            # }else{
            #     ll[i] <- - log(det(omega[,,i])) + sum(diag(cov[,,i]%*%omega[,,i]))
            # }
        }            
        return(ll)
    }
    X0 <- X
    n <- dim(X)[1]
    X <- X[sample(n, n), ]
    ll <- matrix(0, length(v0s), K)
    for(k in 1:K){
        sub <- X[(1:n) %% K != (k-1), ]
        test <-  X[(1:n) %% K == (k-1), ]
        if(copula || sum(is.na(test)) > 0){
            Xtest <- test
        }else{
            Xtest <- NULL
        }
        tmp <- EMGS(sub, v0s, v1, lambda, a, b, epsilon, maxitr, verbose=FALSE, savepath = FALSE, copula, copula_itr, group, Xtest=Xtest, ...)
        if(copula || sum(is.na(test)) > 0){
            cov <- apply(tmp$St, c(1, 2), mean) / dim(test)[1]
            # cov <- tmp$St / dim(test)[1]
        }else{
            cov <- t(test) %*% test / dim(test)[1]
        }
        ll[, k] <- getLL(n = dim(test)[1], cov = cov, omega = tmp$omega, copula) 
    }
    ll.mean <- apply(ll, 1, mean)
    ll.sd <- apply(ll, 1, sd) / sqrt(K - 1)
    v0.min <- v0s[which.min(ll.mean)]
    ll.1se <- ll.mean - (min(ll.mean) + ll.sd[which.min(ll.mean)])
    v0.1se<- v0s[max((1:length(v0s))[ll.1se <= 0])]

    fit.min <- EMGS(X=X0, v0 = v0.min, v1=v1, lambda=lambda, a=a, b=b, epsilon=epsilon, maxitr=maxitr, verbose=FALSE, group = group, copula = copula, copula_itr = copula_itr, ...)
    fit.1se <- EMGS(X=X0, v0 = v0.1se, v1=v1, lambda=lambda, a=a, b=b, epsilon=epsilon, maxitr=maxitr, verbose=FALSE, group = group, copula = copula, copula_itr = copula_itr, ...)
    return(list(ll = ll, v0s = v0s, v0.min = v0.min, v0.1se = v0.1se, fit.min = fit.min, fit.1se = fit.1se))
}


#' Perform Graphical lasso with cross validation
#' 
#' This function performs graphical lasso with cross validation
#' 
#' @param X n by p matrix of data
#' @param rholist list of tunning parameters
#' @param K number of cross validation folds.
#' @param skeptic whether to perform the skeptic transformation
#' @return  fitted objects with largest likelihood and corresponding to 1-SE rule.
#' @examples
#'
cv.glasso <- function(X, rholist = NULL, K = 5, skeptic = FALSE ){
    require("glasso")
    getLL <- function(cov, omega){
        k <- length(omega) 
        ll <- rep(NA, k)
        for(i in 1:k){
            ll[i] <- -log(det(as.matrix(omega[[i]]))) + sum(diag(cov%*%as.matrix(omega[[i]])))
        }
        return(ll)
    }
    n <- dim(X)[1]
    X <- X[sample(n, n), ]
    ll <- matrix(0, length(rholist), K)
    for(k in 1:K){
        sub <- X[(1:n) %% K != (k-1), ]
        test <-  X[(1:n) %% K == (k-1), ]
        if(skeptic){
            sub <- huge.npn(X[(1:n)%%K != (k - 1), ], npn.func = "skeptic")
            test <- huge.npn(X[(1:n)%%K == (k - 1), ], npn.func = "skeptic")
            tmp <- huge.glasso(sub, lambda = rholist)
            score <- getLL(cov = test, omega = tmp$icov)
         }else{
            tmp <- huge.glasso(cov(sub), lambda = rholist)
            cov <- t(test)%*%test/dim(test)[1]
            score <- getLL(cov = cov, omega = tmp$icov)
        }
        ll[,k] <- score
    }
    ll.mean <- apply(ll, 1, mean)
    ll.sd <- apply(ll, 1, sd) / sqrt(K - 1)
    rho.min <- rholist[which.min(ll.mean)]
    ll.1se <- ll.mean - (min(ll.mean) + ll.sd[which.min(ll.mean)])
    rho.1se<- rholist[max((1:length(rholist))[ll.1se <= 0])]

    fit.min <-  huge.glasso(X, lambda = rho.min) 
    fit.1se <- huge.glasso(X, lambda = rho.1se)
    return(list(ll = ll, rholist = rholist, rho.min = rho.min, rho.1se = rho.1se, fit.min = fit.min, fit.1se = fit.1se))
}

 