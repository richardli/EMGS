# srun --pty --partition=short --time=3:00:00 --mem-per-cpu=2500 /bin/bash
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

library(EMGS)
library(tmvtnorm)
library(mvtnorm)
library(glasso)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
library(Rcpp)
library(huge)
source("../codes/BDgraph_sim_func.R")

ns <- c(100, 200, 500)
ps <- c(50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:100
config <- expand.grid(ns, ps, misses, reps)
dim(config)
index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
index <-  which(c(1:dim(config)[1]) %% 200 == index-1)
for(iii in index){
	n <- config[iii, 1]
	p <- config[iii, 2]
	graph <- config[iii, 3]
	itr <- config[iii, 4]
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
                Prec <- fit$Prec
                G <- fit$G
                fit.emgs <- fit$cv.out$fit.min
                fit.glasso <- fit$cv.out.npn$fit.min
                fit.glasso.ric <- fit$out.ric
                fit.glasso.star <- fit$out.star
                fit.glasso.xz <- NULL

                prec.emgs <- fit.emgs$omega[,,1]
                prec.glasso <- as.matrix(fit.glasso$icov[[1]])
                prec.glasso.ric <- as.matrix(fit.glasso.ric$opt.icov)
                prec.glasso.star <- as.matrix(fit.glasso.star$opt.icov)
                prec.glasso.xz <- NULL

                if(!is.null(prec.emgs)){
                    norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star), Prec)
                    roc.emgs <- getROCthres(prec.emgs, G)
                    roc.glasso <- getROCthres(prec.glasso, G)
                    roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                    roc.glasso.star <- getROCthres(prec.glasso.star, G)
                    auc.emgs <- roc.emgs$AUC
                    auc.glasso <- roc.glasso$AUC
                    auc.glasso.ric <- roc.glasso.ric$AUC
                    auc.glasso.star <- roc.glasso.star$AUC
                    graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                    graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                    graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                    graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                    f1 <- getF1(list(graph.emgs, graph.glasso, graph.glasso.ric, graph.glasso.star), G)

                }else{
                    norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star), Prec)
                    for(ii in 1:length(norms)) norms[[ii]] <- rbind(NA, norms[[ii]])
                    roc.glasso <- getROCthres(prec.glasso, G)
                    roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                    roc.glasso.star <- getROCthres(prec.glasso.star, G)
                    auc.emgs <- NA
                    auc.glasso <- roc.glasso$AUC
                    auc.glasso.ric <- roc.glasso.ric$AUC
                    auc.glasso.star <- roc.glasso.star$AUC
                    graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                    graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                    graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                    graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                    f1 <- c(NA, getF1(list(graph.glasso, graph.glasso.ric, graph.glasso.star), G))
                }

                # AUC
                AUC_lam <- c(getROCpath(fit$out$path, G)$AUC, getROCpath(fit$out.npn$path, G)$AUC, NA, NA)
                # F1
                graph.raw <- fit$cv.out$fit.min$path[,,1]
                graph.gl  <- prec.glasso
                graph.gl[graph.gl != 0] <- 1        
                graph.ric  <- prec.glasso.ric
                graph.ric[graph.ric != 0] <- 1 
                graph.star <- prec.glasso.star
                graph.star[graph.star != 0] <- 1 
                F1_raw <-  getF1(list(graph.raw, graph.gl, graph.ric, graph.star), G)

                eval.g <- list(norms=norms, auc = c(auc.emgs, auc.glasso, auc.glasso.ric, auc.glasso.star), f1 = f1, Auc_lam = AUC_lam, f1_raw = F1_raw)
	}else{
        eval.g <- NULL
    }
	if(file.exists(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed-1.rda", sep="_"))){
        load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed-1.rda", sep="_"))
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
            save(fit, file = paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed-1.rda", sep="_"))
        }
                    Prec <- fit$Prec
                    G <- fit$G
                    fit.emgs <- fit$cv.out$fit.min
                    fit.glasso <- fit$cv.out.npn$fit.min
                    fit.glasso.ric <- fit$out.ric
                    fit.glasso.star <- fit$out.star
                    fit.glasso.xz <- fit$cv.out.xz$fit.min

                    prec.emgs <- fit.emgs$omega[,,1]
                    prec.glasso <- as.matrix(fit.glasso$icov[[1]])
                    prec.glasso.ric <- as.matrix(fit.glasso.ric$opt.icov)
                    prec.glasso.star <- as.matrix(fit.glasso.star$opt.icov)
                    prec.glasso.xz <- as.matrix(fit.glasso.xz$icov[[1]])

                    if(!is.null(prec.emgs)){
                        norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star, prec.glasso.xz), Prec)
                        roc.emgs <- getROCthres(prec.emgs, G)
                        roc.glasso <- getROCthres(prec.glasso, G)
                        roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                        roc.glasso.star <- getROCthres(prec.glasso.star, G)
                        roc.glasso.xz <- getROCthres(prec.glasso.xz, G)
                        auc.emgs <- roc.emgs$AUC
                        auc.glasso <- roc.glasso$AUC
                        auc.glasso.ric <- roc.glasso.ric$AUC
                        auc.glasso.star <- roc.glasso.star$AUC
                        auc.glasso.xz <- roc.glasso.xz$AUC
                        graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                        graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                        graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                        graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                        graph.glasso.xz <- Threshold(prec.glasso.xz, sum(G)/2)
                        f1 <- getF1(list(graph.emgs, graph.glasso, graph.glasso.ric, graph.glasso.star, graph.glasso.xz), G)
                        # AUC
                        AUC_lam <- c(getROCpath(fit$out$path, G)$AUC, getROCpath(fit$out.npn$path, G)$AUC, NA, NA, getROCpath(fit$out.xz$path, G)$AUC)
                        # F1
                        graph.raw <- fit$cv.out$fit.min$path[,,1]
                        graph.gl  <- prec.glasso
                        graph.gl[graph.gl != 0] <- 1        
                        graph.ric  <- prec.glasso.ric
                        graph.ric[graph.ric != 0] <- 1 
                        graph.star <- prec.glasso.star
                        graph.star[graph.star != 0] <- 1 
                        graph.xz <- prec.glasso.xz
                        graph.xz[graph.xz != 0] <- 1 
                        F1_raw <-  getF1(list(graph.raw, graph.gl, graph.ric, graph.star, graph.xz), G)

                    }else{
                        norms <- getNorms(list(prec.glasso, prec.glasso.ric, prec.glasso.star, prec.glasso.xz), Prec)
                        for(ii in 1:length(norms)) norms[[ii]] <- rbind(NA, norms[[ii]])
                        # roc.emgs <- getROCthres(prec.emgs, G)
                        roc.glasso <- getROCthres(prec.glasso, G)
                        roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                        roc.glasso.star <- getROCthres(prec.glasso.star, G)
                        roc.glasso.xz <- getROCthres(prec.glasso.xz, G)
                        auc.emgs <- NA
                        auc.glasso <- roc.glasso$AUC
                        auc.glasso.ric <- roc.glasso.ric$AUC
                        auc.glasso.star <- roc.glasso.star$AUC
                        auc.glasso.xz <- roc.glasso.xz$AUC
                        # graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                        graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                        graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                        graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                        graph.glasso.xz <- Threshold(prec.glasso.xz, sum(G)/2)
                        f1 <- c(NA, getF1(list(graph.glasso, graph.glasso.ric, graph.glasso.star, graph.glasso.xz), G))
                        # AUC
                        AUC_lam <- c(NA, getROCpath(fit$out.npn$path, G)$AUC, NA, NA, getROCpath(fit$out.xz$path, G)$AUC)
                        # F1
                        graph.gl  <- prec.glasso
                        graph.gl[graph.gl != 0] <- 1        
                        graph.ric  <- prec.glasso.ric
                        graph.ric[graph.ric != 0] <- 1 
                        graph.star <- prec.glasso.star
                        graph.star[graph.star != 0] <- 1 
                        graph.xz <- prec.glasso.xz
                        graph.xz[graph.xz != 0] <- 1 
                        F1_raw <-  c(NA, getF1(list(graph.gl, graph.ric, graph.star, graph.xz), G))     
                    }

                    eval.m <- list(norms=norms, auc = c(auc.emgs, auc.glasso, auc.glasso.ric, auc.glasso.star, auc.glasso.xz), f1 = f1, Auc_lam = AUC_lam, f1_raw = F1_raw)
    }else{
        eval.m <- NULL
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
            Prec <- fit$Prec
            G <- fit$G
            fit.emgs <- fit$cv.out$fit.min
            fit.glasso <- fit$cv.out.npn$fit.min
            fit.glasso.ric <- fit$out.ric
            fit.glasso.star <- fit$out.star
            fit.glasso.xz <- fit$cv.out.xz$fit.min

            prec.emgs <- fit.emgs$omega[,,1]
            prec.glasso <- as.matrix(fit.glasso$icov[[1]])
            prec.glasso.ric <- as.matrix(fit.glasso.ric$opt.icov)
            prec.glasso.star <- as.matrix(fit.glasso.star$opt.icov)
            prec.glasso.xz <- as.matrix(fit.glasso.xz$icov[[1]])

            if(!is.null(prec.emgs)){
                norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star, prec.glasso.xz), Prec)
                roc.emgs <- getROCthres(prec.emgs, G)
                roc.glasso <- getROCthres(prec.glasso, G)
                roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                roc.glasso.star <- getROCthres(prec.glasso.star, G)
                roc.glasso.xz <- getROCthres(prec.glasso.xz, G)
                auc.emgs <- roc.emgs$AUC
                auc.glasso <- roc.glasso$AUC
                auc.glasso.ric <- roc.glasso.ric$AUC
                auc.glasso.star <- roc.glasso.star$AUC
                auc.glasso.xz <- roc.glasso.xz$AUC
                graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                graph.glasso.xz <- Threshold(prec.glasso.xz, sum(G)/2)
                f1 <- getF1(list(graph.emgs, graph.glasso, graph.glasso.ric, graph.glasso.star, graph.glasso.xz), G)
                # AUC
                AUC_lam <- c(getROCpath(fit$out$path, G)$AUC, getROCpath(fit$out.npn$path, G)$AUC, NA, NA, getROCpath(fit$out.xz$path, G)$AUC)
                # F1
                graph.raw <- fit$cv.out$fit.min$path[,,1]
                graph.gl  <- prec.glasso
                graph.gl[graph.gl != 0] <- 1        
                graph.ric  <- prec.glasso.ric
                graph.ric[graph.ric != 0] <- 1 
                graph.star <- prec.glasso.star
                graph.star[graph.star != 0] <- 1 
                graph.xz <- prec.glasso.xz
                graph.xz[graph.xz != 0] <- 1 
                F1_raw <-  getF1(list(graph.raw, graph.gl, graph.ric, graph.star, graph.xz), G)

            }else{
                norms <- getNorms(list(prec.glasso, prec.glasso.ric, prec.glasso.star, prec.glasso.xz), Prec)
                for(ii in 1:length(norms)) norms[[ii]] <- rbind(NA, norms[[ii]])
                # roc.emgs <- getROCthres(prec.emgs, G)
                roc.glasso <- getROCthres(prec.glasso, G)
                roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
                roc.glasso.star <- getROCthres(prec.glasso.star, G)
                roc.glasso.xz <- getROCthres(prec.glasso.xz, G)
                auc.emgs <- NA
                auc.glasso <- roc.glasso$AUC
                auc.glasso.ric <- roc.glasso.ric$AUC
                auc.glasso.star <- roc.glasso.star$AUC
                auc.glasso.xz <- roc.glasso.xz$AUC
                # graph.emgs <- Threshold(prec.emgs, sum(G)/2)
                graph.glasso <- Threshold(prec.glasso, sum(G)/2)
                graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
                graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
                graph.glasso.xz <- Threshold(prec.glasso.xz, sum(G)/2)
                f1 <- c(NA, getF1(list(graph.glasso, graph.glasso.ric, graph.glasso.star, graph.glasso.xz), G))
                # AUC
                AUC_lam <- c(NA, getROCpath(fit$out.npn$path, G)$AUC, NA, NA, getROCpath(fit$out.xz$path, G)$AUC)
                # F1
                graph.gl  <- prec.glasso
                graph.gl[graph.gl != 0] <- 1        
                graph.ric  <- prec.glasso.ric
                graph.ric[graph.ric != 0] <- 1 
                graph.star <- prec.glasso.star
                graph.star[graph.star != 0] <- 1 
                graph.xz <- prec.glasso.xz
                graph.xz[graph.xz != 0] <- 1 
                F1_raw <-  c(NA, getF1(list(graph.gl, graph.ric, graph.star, graph.xz), G))     
            }

            eval.m <- list(norms=norms, auc = c(auc.emgs, auc.glasso, auc.glasso.ric, auc.glasso.star, auc.glasso.xz), f1 = f1, Auc_lam = AUC_lam, f1_raw = F1_raw)
    }else{
        eval.m <- NULL
    }
    out <- list(Gaussian = eval.g, Mixed = eval.m)
    save(out, file = paste("../data/metrics/", n, p, graph, itr, "metrics.rda", sep="_"))     
	cat("=======================FIN=======================")
}
