# Help function to simulate various types of graphs
BDprec <- function(P, type){
	wi <- diag(P)
	g <- matrix(0, P, P)
	if(type == "circle"){
		for(i in 2:P){
			wi[i, i-1] <- wi[i-1, i] <- 0.5
			g[i, i-1] <- g[i-1, i] <- 1
		}
		wi[1, P] <- wi[P, 1] <- 0.4
		g[1, P] <- g[P, 1] <- 1
	}else if(type == "star"){
		for(i in 2:P){
			wi[1, i] <- wi[i, 1] <- 0.1
			g[1, i] <- g[i, 1] <- 1
		}
	}else if(type == "AR1"){
		# default to be sigma_ij = 0.7^|i-j|
		tmp <- bdgraph.sim(n = 1, p = P, vis = FALSE, graph = "AR1" )
		wi <- tmp$K
		g <- tmp$G
	}else if(type == "AR2"){
		# default to be 1, .5, .25
		tmp <- bdgraph.sim(n = 1, p = P, vis = FALSE, graph = "AR2" )
		wi <- tmp$K
		g <- tmp$G
	}else if(type == "random"){
		tmp <- bdgraph.sim(n = 1, p = P, vis = FALSE, graph = "random", prob = 2/(P-1))
		wi <- tmp$K
		g <- tmp$G
	}else if(type == "cluster"){
		tmp <- bdgraph.sim(n = 1, p = P, vis = FALSE, graph = "cluster")
		wi <- tmp$K
		g <- tmp$G
	}else if(type == "scale-free"){
		tmp <- bdgraph.sim(n = 1, p = P, vis = FALSE, graph = "scale-free")
		wi <- tmp$K
		g <- tmp$G
	}
	return(list(wi = wi, g = g))
}

# Help function to calculate the posterior conditional inclusion probability
getp <- function(v0, v1, prec, tau=1){
	expit <- function(x){
		p1 <- exp(x) / (1 + exp(x))
		p2 <- 1 / (1 + exp(-x))
		p1[is.na(p1)] <- 0
		p2[is.na(p2)] <- 0
		return(p1 * (x < 0) + p2 * (x >= 0))	
	}
	a <- v1 / v0
	prob <- expit( prec^2 / (2 * v0^2/tau * a^2) * (a^2 - 1) - log(a))
	return(prob)
}

# Help function to calculate F1-score
getF1 <- function(paths, truth){
	f1 <- rep(0, length(paths))
	truth[truth != 0] <- 1
	diag(truth) <- NA
	truth[upper.tri(truth)] <- NA
	for(i in 1:length(paths)){
		pred <- as.matrix(paths[[i]])
		pred[pred != 0] <- 1
		diag(pred) <- NA
		pred[upper.tri(pred)] <- NA
		TP <- length(intersect(which(pred == 1), which(truth == 1)))
		FP <- length(intersect(which(pred == 1), which(truth == 0)))
		FN <- length(intersect(which(pred == 0), which(truth == 1)))
		f1[i] <- 2 * TP / (2 * TP + FP + FN)
	}
	return(f1)
}

###########################################################################
get4norm <- function(list, cov){
	mn <- NULL
	sn <- NULL
	fn <- NULL
	infn <- NULL
	for(i in 1:length(list)){
		mtemp <- norm(as.matrix(list[[i]]) - cov, type = "m")
		stemp <- base::norm(as.matrix(list[[i]]) - cov, type = "2")
		ftemp <- norm(as.matrix(list[[i]]) - cov, type = "f")
		inftemp <- max(apply(as.matrix(list[[i]]) - cov, 1, function(x){sum(abs(x))}))
		mn <- c(mn, mtemp)
		sn <- c(sn, stemp)
		fn <- c(fn, ftemp)
		infn <- c(infn, inftemp)
	}
	out <- cbind(mn, sn, fn, infn)
	colnames(out) <- c("M norm", "S norm", "F norm", "inf norm")
	return(out)
}
###########################################################################
pathplot <- function(v0s, path, G, normalize = FALSE, xlab="", ylab = "", main = ""){
	if(normalize){
		for(i in 1:dim(path)[1]){ path[i, , ] <- cov2cor(path[i, , ]) }
	}
	for(i in 1:dim(path)[1]){ 
		diag(path[i, , ]) <- NA 
	}
	plot(v0s, rep(0, length(v0s)), ylim = range(path, na.rm = TRUE), col = "white", xlab = xlab, main = main, ylab = ylab)
	M <- dim(G)[1]
	for(i in 2:M){
		for(j in 1:(i-1)){
			if(G[i, j] != 0){
				lines(v0s, path[, i, j], type = "l", col = "darkblue")
			}
			if(G[i, j] == 0){
				lines(v0s, path[, i, j], type = "l", col = "darkgreen", lty = 4)
			}
		}
	}
	legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("True discovery", "False discovery"))
}
###########################################################################
getROC <- function(path, G){
	for(i in 1:dim(path)[1]){ 
		diag(path[i, , ]) <- NA 
	}
	diag(G) <- NA
	tfpr <- matrix(0, dim(path)[1], 2)
	colnames(tfpr) <- c("TPR", "FPR")
	for(i in 1:dim(path)[1]){
		pred <-  path[i, , ] != 0
		tfpr[i, 1] <- length(intersect(which(pred == 1), which(G == 1))) / length(which(G == 1))
		tfpr[i, 2] <- length(intersect(which(pred == 1), which(G == 0))) / length(which(G == 0))
	}
	tfpr <- rbind(c(0, 0), tfpr, c(1, 1))
	tfpr <- tfpr[order(tfpr[, 2], tfpr[, 1]), ]
	return(tfpr)
}

###########################################################################
Threshold <- function(prec, nedge){
	tmp <- prec
	diag <- diag(prec)
	tmp[!lower.tri(tmp)] <- NA
	thre <- sort(abs(as.numeric(tmp)), decreasing = TRUE)[nedge]
	prec[abs(prec) < thre] <- NA
	diag(prec) <- diag
	return(prec)
}

###########################################################################
getROCthres <- function(prec, G){
	diag(prec) <- NA
	diag(G) <- NA
	# cut <- sort(abs(prec))
	# tfpr <- matrix(0, length(cut), 2)
	# colnames(tfpr) <- c("TPR", "FPR")
	# for(i in 1:length(cut)){
	# 	pred <-  prec
	# 	pred[abs(prec) <= cut[i]] <- 0
	# 	tfpr[i, 1] <- length(intersect(which(pred != 0), which(G != 0))) / length(which(G == 1))
	# 	tfpr[i, 2] <- length(intersect(which(pred != 0), which(G == 0))) / length(which(G == 0))
	# }
	# tfpr <- rbind(c(0, 0), tfpr, c(1, 1))
	# tfpr <- tfpr[order(tfpr[, 2], tfpr[, 1]), ]
	# return(tfpr)
	pred <- prediction(as.numeric(abs(prec)), as.numeric(G))
	tfpr <- cbind(TPR = performance(pred, "tpr")@y.values[[1]],
					FPR = performance(pred, "fpr")@y.values[[1]]) 
	AUC <- performance(pred, "auc")@y.values[[1]]
	return(list(tfpr=tfpr, AUC=AUC))
}

###########################################################################
getROCpath <- function(path, G){
	if(length(dim(path)) == 3){
		tmp <- vector("list", dim(path)[3])
		for(i in 1:length(tmp)) tmp[[i]] <- path[,,i]
	}else{
		tmp <- path
	}
	tmp[[length(tmp) + 1]] <- matrix(0, dim(tmp[[1]])[1], dim(tmp[[1]])[2])
	tmp[[length(tmp) + 1]] <- matrix(1, dim(tmp[[1]])[1], dim(tmp[[1]])[2])
	roc <- huge.roc(tmp, theta = G, verbose = FALSE)
	return(roc)
}
