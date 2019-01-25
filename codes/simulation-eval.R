library(EMGS)
library(tmvtnorm)
library(mvtnorm)
library(glasso)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
library(Rcpp)
source("../codes/BDgraph_sim_func.R")
ns <- c(100, 200, 500)
ps <- c(50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:100
norms <- c("M norm", "S norm", "F norm", "Inf norm", "M norm0", "S norm0", "F norm0", "Inf norm0", "AUC", "F1" , "AUC-lam", "F1-raw" )
mets <- array(NA, dim = c(5, length(norms), 3, 3, 4, 100))
dimnames(mets)[[1]] <- c("EMGS", "Glasso-thres", "Glasso-ric", "Glasso-star", "Glasso-XZ")
dimnames(mets)[[2]] <- norms
dimnames(mets)[[3]] <- paste("n=", ns)
dimnames(mets)[[4]] <- paste("p=", ps)
dimnames(mets)[[5]] <- misses
mets.m <- mets

for(i in 1:3){
for(j in 1:3){
for(k in 1:4){
for(l in 1:100){
	n <- ns[i]
	p <- ps[j]
	graph <- misses[k]
	itr <- reps[l]
	cal.g <- cal.m <- TRUE
	if(file.exists(paste("../data/metrics/", n, p, graph, itr, "metrics.rda", sep="_"))){
		load(paste("../data/metrics/", n, p, graph, itr, "metrics.rda", sep="_"))
		if(!is.null(out$Gaussian) && dim(out$Gaussian$norms$norm.prec)[1] > 2){
			cal.g <- FALSE
			eval.g <- out$Gaussian
			mets[1:4, , i,j,k,l] <- cbind(eval.g$norms$norm.prec, eval.g$norms$norm.prec, eval.g$auc, eval.g$f1, eval.g$Auc_lam, eval.g$f1_raw)
		}
		if(!is.null(out$Mixed) && dim(out$Mixed$norms$norm.prec)[1] > 2){
			eval.m <- out$Mixed
			cal.m <- FALSE
			mets.m[, , i,j,k,l] <- cbind(eval.m$norms$norm.prec, eval.m$norms$norm.prec, eval.m$auc, eval.m$f1, eval.m$Auc_lam, eval.m$f1_raw)
		}
	} 

	if(cal.g && file.exists(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))){
		###################################################################
		## Gaussian case
		load(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
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
		mets[1:4, , i,j,k,l] <- cbind(eval.g$norms$norm.prec, eval.g$norms$norm.prec, eval.g$auc, eval.g$f1, AUC_lam, F1_raw)
	}else{
		eval.g <- NULL
	}

	###################################################################
	## Mixed case
	if(cal.m && file.exists(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))){
		load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
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
		mets.m[, , i,j,k,l] <- cbind(eval.m$norms$norm.prec, eval.m$norms$norm.prec, eval.m$auc, eval.m$f1, AUC_lam, F1_raw)
	}else{
		eval.m <- NULL
	}

	###################################################################
	out <- list(Gaussian = eval.g, Mixed = eval.m)
	save(out, file = paste("../data/metrics/", n, p, graph, itr, "metrics.rda", sep="_"))
	}
	print(paste("Gaussian", n, p, graph))
	print(apply(mets[, , i, j, k, ], c(1, 2), mean, na.rm = TRUE))
	print(paste("Mixed", n, p, graph))
	print(apply(mets.m[, , i, j, k, ], c(1, 2), mean, na.rm = TRUE))

}}}



m0 <- apply(mets, c(1, 2, 3, 4, 5), mean, na.rm = TRUE)
dimnames(m0) 



library(xml2)
library(knitr)
library(kableExtra)
library(plyr)
norms <- c( "S norm", "F norm", "AUC", "F1" , "AUC-lam", "F1-raw" )
config.c <- expand.grid( ns, norms)
config.r <- expand.grid( dimnames(mets)[[1]], ps[2], misses)

tab <- matrix(NA, dim(config.r)[1], dim(config.c)[1]+3)
tab <- data.frame(tab)
tab[, c(1, 2, 3)] <- ""
for(i in 1:dim(tab)[1]){
	for(j in 1:(dim(tab)[2]-3)){
		ii = which(dimnames(m0)[[1]] == config.r[i,1])
		jj = which(dimnames(m0)[[2]] == config.c[j,2])
		kk = which(dimnames(m0)[[3]] == paste("n=", config.c[j, 1]))
		ll = which(dimnames(m0)[[4]] == paste("p=", config.r[i, 2]))
		uu = which(dimnames(m0)[[5]] == config.r[i,3])
		tab[i, j+3] <- m0[ii, jj, kk, ll, uu]
	}
}
config.r0 <- config.r
for(i in 2:dim(config.r)[1]){
	if(config.r[i, 2] == config.r[i-1, 2]) config.r0[i, 2] <- NA
}
for(i in 2:dim(config.r)[1]){
	if(config.r[i, 3] == config.r[i-1, 3]) config.r0[i, 3] <- NA
}
tab[,1:3] <- config.r0[, c(3, 2, 1)]
tab <- tab[config.r[, 1] %in% c("Glasso-XZ") == FALSE, ]
colnames(tab) <- c("graph", "p", "method", config.c[,1])
tab[,3] <- revalue(tab[,3], c("Glasso-thres" = "GL", "Glasso-ric" = "GL-RIC", "Glasso-star"="GL-StARS"))

options(knitr.kable.NA = '')
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the precision matrix for Graphical model estimation.", align = "c", row.names = FALSE, linesep = "") %>% add_header_above(c(" ", " ", " ", "S"=3, "F"=3, "AUC"=3, "F1" = 3, "AUC"=3, "F1" = 3)) %>% kable_styling(latex_options = "hold_position", font_size = 7)

