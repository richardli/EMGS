cv.glasso.xz <- function (X, rholist = NULL, K = 5){
    require("glasso")
    getLL <- function(Cov, omega) {
        k <- length(omega)
        ll <- rep(NA, k)
        for (i in 1:k) {
            ll[i] <- -log(det(as.matrix(omega[[i]]))) + sum(diag(Cov %*% as.matrix(omega[[i]])))
        }
        return(ll)
    }
    if (is.null(rholist)) {
        lambda.min.ratio = 0.1
        d <- dim(X)[2]
        lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
        lambda.min = lambda.min.ratio * lambda.max
        rholist = rev(exp(seq(log(lambda.max), log(lambda.min), 
            length = 30)))
    }
    n <- dim(X)[1]
    X <- X[sample(n, n), ]
    ll <- matrix(0, length(rholist), K)
    for (k in 1:K) {
        sub <- huge.npn(X[(1:n)%%K != (k - 1), ], npn.func = "skeptic")
        test <- huge.npn(X[(1:n)%%K == (k - 1), ], npn.func = "skeptic")
        tmp <- huge.glasso(sub, lambda = rholist)
        score <- getLL(Cov = test, omega = tmp$icov)
        ll[, k] <- score
    }
    ll.mean <- apply(ll, 1, mean)
    ll.sd <- apply(ll, 1, sd)/sqrt(K - 1)
    rho.min <- rholist[which.min(ll.mean)]
    ll.1se <- ll.mean - (min(ll.mean) + ll.sd[which.min(ll.mean)])
    rho.1se <- rholist[max((1:length(rholist))[ll.1se <= 0])]
    fit.min <- huge.glasso(X, lambda = rho.min)
    fit.1se <- huge.glasso(X, lambda = rho.1se)
    return(list(ll = ll, rholist = rholist, rho.min = rho.min, 
        rho.1se = rho.1se, fit.min = fit.min, fit.1se = fit.1se))
}



ns <- c(100, 200, 500)
ps <- c(50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:100


for(i in 1:3){
for(j in 1:3){
for(k in 1:4){
for(l in 1:100){
	n <- ns[i]
	p <- ps[j]
	graph <- misses[k]
	itr <- reps[l]
	if(file.exists(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))){
		load(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
		if(is.null(fit$out.ric)){
			load(paste("../data/simG/graphsim", n, p, graph, itr, "gaussian.rda", sep="_"))
			X <- sim$data
			out1 = huge.select(fit$out.npn, criterion = "ric")
			out2 = huge.select(fit$out.npn, criterion = "stars")
			fit$out.ric <- out1
			fit$out.star <- out2
			save(fit, file = paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))			
		}
	}
	if(file.exists(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))){
		load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
		if(is.null(fit$out.ric)){
			load(paste("../data/simM/graphsim", n, p, graph, itr, "mixed.rda", sep="_"))
			X <- sim$data
			out1 = huge.select(fit$out.npn, criterion = "ric")
			out2 = huge.select(fit$out.npn, criterion = "stars")
			XX <- huge.npn(X, npn.func = "skeptic")
			out3 <- huge(XX, method = "glasso", nlambda=40)
			cv.out3 = cv.glasso.xz(X, rholist = out3$lambda, K=5)

			fit$out.ric <- out1
			fit$out.star <- out2
			fit$out.xz <- out3
			fit$cv.out.xz <- cv.out3
			save(fit, file = paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
		}
	}
	cat(".")
}
cat("\n")
cat(paste(n, p, graph, "done\n"))
}
}
}