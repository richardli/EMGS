#' Calculate the matrix norms
#' 
#' This function calculates the M/S/F/Inf norms for the bias of a list of estimated precision matrices with the truth. 
#' 
#' @param obj a fitted EMGS object or a list of precision matrices
#' @param prec true precision matrix
#' @param cov true covariance matrix
#' @return norms of the bias for (1) precision matrix, (2) standardized precision matrix, (3) covariance matrix, and (4) correlation matrix
#' @examples
#' 


getNorms <- function(obj, prec, cov = NULL){
	
	getone <- function(est, truth){
		mtemp <- norm(as.matrix(est) - truth, type = "m")
		stemp <- base::norm(as.matrix(est) - truth, type = "2")
		ftemp <- norm(as.matrix(est) - truth, type = "f")
		inftemp <- max(apply(as.matrix(est) - truth, 1, function(x){sum(abs(x))}))
		return(c(mtemp, stemp, ftemp, inftemp))
	}

	if(is.null(cov)) cov = solve(prec)
	if(class(obj) == "list"){
		out1 <- out2 <- out3 <- out4 <- matrix(NA, length(obj), 4)
		for(i in 1:length(obj)){
			out1[i, ] <- getone(obj[[i]], prec);
			out2[i, ] <- getone(cov2cor(obj[[i]]), cov2cor(prec));
			out3[i, ] <- getone(solve(obj[[i]]), cov);
			out4[i, ] <- getone(cov2cor(solve(obj[[i]])), cov2cor(cov));
		}
	}else if(class(obj) == "emgs"){
		# fitted object
		out1 <- out2 <- out3 <- out4 <- matrix(NA, dim(obj$omega)[3], 4)
		for(i in 1:dim(obj$omega)[3]){
			out1[i, ] <- getone(obj$omega[, , i], prec);
			out2[i, ] <- getone(cov2cor(obj$omega[, , i]), cov2cor(prec));
			out3[i, ] <- getone(solve(obj$omega[, , i]), cov);
			out4[i, ] <- getone(cov2cor(solve(obj$omega[, , i])), cov2cor(cov));
		}
	}else if(class(obj) == "array"){
		out1 <- out2 <- out3 <- out4 <- matrix(NA, dim(obj)[3], 4)
		for(i in 1:dim(obj)[3]){
			out1[i, ] <- getone(obj[, , i], prec);
			out2[i, ] <- getone(cov2cor(obj[, , i]), cov2cor(prec));
			out3[i, ] <- getone(solve(obj[, , i]), cov);
			out4[i, ] <- getone(cov2cor(solve(obj[, , i])), cov2cor(cov));
		}
	}
	colnames(out1) <- colnames(out2) <- colnames(out3) <- colnames(out4) <- c("M norm", "S norm", "F norm", "inf norm")

	out <- list(norm.prec = out1, 
				norm.prec.std = out2, 
				norm.cov = out3, 
				norm.corr = out4)

	return(out)
}

#' Plot the regularization path
#' 
#' This function plots the regularization path
#' 
#' @param v0s the vector of parameter v0
#' @param obj fitted EMGS object, or a three-dimensional array of precision matrices
#' @param thres threshold of slab
#' @param normalize logical indicator to plot precision matrices or partial correlation matrices
#' @param xlab x axis label
#' @param ylab y axis label
#' @param main plot title
#' @return the regularization path
#' @examples
#' 


path.plot <- function(v0s, obj, G, thres = 0.5, normalize = FALSE, xlab="", ylab = "", main = ""){
	if(class(obj) == "emgs"){
		prob <- obj$EZ
		obj <- obj$omega
	}else{
		prob <- NULL
	}
	if(normalize){
		for(i in 1:dim(obj)[3]){ 
			obj[, , i] <- -cov2cor(obj[, , i]) 
		}
	}
	if(!is.null(prob)){
		for(i in 1:dim(obj)[3]){ 
			obj[, , i] <- obj[, , i] * (prob[, , i] > thres)
		}
	}
	for(i in 1:dim(obj)[3]){ 
		diag(obj[, , i]) <- NA 
	}
	plot(v0s, rep(0, length(v0s)), ylim = range(obj, na.rm = TRUE), col = "white", xlab = xlab, main = main, ylab = ylab)
	M <- dim(G)[1]
	for(i in 2:M){
		for(j in 1:(i-1)){
			if(G[i, j] == 0){
				lines(v0s, obj[i, j, ], type = "l", col = "darkgreen", lty = 4)
			}
		}
	}
	for(i in 2:M){
		for(j in 1:(i-1)){
			if(G[i, j] != 0){
				lines(v0s, obj[i, j, ], type = "l", col = "darkblue")
			}
		}
	}
	legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("True positive", "False positive"))
}

#' Calculate the true positive rate and false positive rate
#' 
#' This function calculates the true positive rate and false positive rate of edge discovery
#' 
#' @param obj fitted EMGS object, or a three-dimensional array of precision matrices
#' @param G the true graph. Diagonal entries can take any value.
#' @param thres threshold of slab
#' @return a matrix of TPR and FPR
#' @examples
#' 

get.ROC <- function(obj, G, thres = 0.5){
	if(class(obj) == "emgs"){
		obj <- obj$EZ > thres
	}
	for(i in 1:dim(obj)[3]){ 
		diag(obj[, , i]) <- NA 
	}
	diag(G) <- NA
	tfpr <- matrix(0, dim(obj)[3], 2)
	colnames(tfpr) <- c("TPR", "FPR")
	for(i in 1:dim(obj)[3]){
		pred <-  obj[, , i] != 0
		tfpr[i, 1] <- length(intersect(which(pred == 1), which(G != 0))) / max(length(which(G != 0)), 1)
		tfpr[i, 2] <- length(intersect(which(pred == 1), which(G == 0))) / max(length(which(G == 0)), 1)
	}
	tfpr <- rbind(c(0, 0), tfpr, c(1, 1))
	tfpr <- tfpr[order(tfpr[, 2], tfpr[, 1]), ]
	return(tfpr)
}

#' Calculate the true positive rate and false positive rate
#' 
#' This function calculates the true positive rate and false positive rate of edge discovery
#' 
#' @param obj fitted EMGS object, or a three-dimensional array of precision matrices
#' @param G the true graph. Diagonal entries can take any value.
#' @param thres threshold of slab
#' @return a matrix of TPR and FPR
#' @examples
#' 
get.ROC.thre <- function(obj, G, verbose = FALSE){
	gcinfo(verbose = FALSE)
    ROC = list()
    G = as.matrix(G)
    d = ncol(G)
    pos.total = sum(G != 0)
    neg.total = d * (d - 1) - pos.total
    if (verbose) 
        cat("Adapted: Computing F1 scores, false positive rates and true positive rates....")
    K <- length(unique(as.vector(obj)))
    cutoff <- sort(unique(as.vector(obj)), decreasing = TRUE)
    ROC$tp = rep(0, K)
    ROC$fp = rep(0, K)
    ROC$F1 = rep(0, K)
    for (r in 1:K) {
        tmp = as.matrix(obj > cutoff[r])
        tp.all = (G != 0) * (tmp != 0)
        diag(tp.all) = 0
        ROC$tp[r] <- sum(tp.all != 0)/pos.total
        fp.all = (G == 0) * (tmp != 0)
        diag(fp.all) = 0
        ROC$fp[r] <- sum(fp.all != 0)/neg.total
        fn = 1 - ROC$tp[r]
        precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
        recall = ROC$tp[r]/(ROC$tp[r] + fn)
        ROC$F1[r] = 2 * precision * recall/(precision + recall)
        if (is.na(ROC$F1[r])) 
            ROC$F1[r] = 0
    }
    if (verbose) 
        cat("done.\n")
    rm(precision, recall, tp.all, fp.all, obj, G, fn)
    gc()
    ord.fp = order(ROC$fp)
    tmp1 = ROC$fp[ord.fp]
    tmp2 = ROC$tp[ord.fp]
    # par(mfrow = c(1, 1))
    # plot(tmp1, tmp2, type = "b", main = "ROC Curve", xlab = "False Postive Rate", 
    #     ylab = "True Postive Rate", ylim = c(0, 1))
    ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
    rm(ord.fp, tmp1, tmp2)
    gc()
    class(ROC) = "roc"
    return(ROC)

}
