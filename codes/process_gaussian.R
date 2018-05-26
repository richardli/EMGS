
##########################################################################
## Process results
##########################################################################
remove(list = ls())
library(EMGS)
library(tmvtnorm)
library(mvtnorm)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
library(Rcpp)
source("../codes/BDgraph_sim_func.R")
config <- NULL
ns <- 50 #c(50, 100, 500)
ps <- 50 #c(10, 50)
simtypes <- "random" # c("AR1", "AR2", "random", "cluster")
fittype <- c("gaussian", "mixed")[1]
config <- expand.grid(ns, ps, simtypes, stringsAsFactors = FALSE)
v0s <- seq(0.001, 0.2, len = 40)

files <- list.files("../data/fit")
alllist <- vector("list", dim(config)[1])
for(index in 1:dim(config)[1]){
	args_in <- as.character(config[index, ])
	p = as.numeric(args_in[2])
	n = as.numeric(args_in[1])
	graph = args_in[3]
	itr = as.numeric(args_in[4])
	tmp <- gsub(paste("graphfit", n, p, graph, "", sep="_"), "!", files)
	tmp <- gsub(paste0("_", fittype, ".rda"), "!", tmp)
	tmp <- unlist(strsplit(tmp, "!"))
	# force the non-rep number to be NA
	seeds <- sort(unique(as.numeric(tmp)))
	ave.curve <- ave.curve2 <- rep(0, p^2/2)
	ave.norms <- matrix(0, 2, 4)
	ave.f1 <- rep(0, 2)
	all.norms <- array(NA, c(length(seeds), 2, 4))
	all.f1 <- matrix(NA, length(seeds), 2)
	rownames(ave.norms) <- names(ave.f1) <- dimnames(all.norms)[[2]] <- colnames(all.f1) <- c("EMGS", "Glasso")

	count <- 0
	for(index0 in 1:length(seeds)){
		itr <- seeds[index0]		
		load( paste("../data/fit/graphfit", n, p, graph, itr, paste0(fittype, ".rda"), sep="_"))
		load( paste("../data/sim/graphsim", n, p, graph, itr, paste0(fittype, ".rda"), sep="_"))
		if(sum(sim$G) < 1){next}
		count <- count + 1

		Prec <- sim$K  
		X <- sim$data
		out <- fit$out
		out.npn <- fit$out.npn
		nedges <- apply(out$path, 3, sum)
		mindist <- which.min(abs(nedges - sum(sim$G)))
		# getNorms(list(out$omega[,,mindist]), Prec)[[1]]
		# getF1(list(out$path[,,mindist]), sim$G)
		# path.plot(v0s, out, sim$G, thres = 0.5, normalize=TRUE, main = paste("N =",  dim(X)[1], "M =", dim(X)[2]))

		nedges.huge <- rep(0, length(out.npn$path))
		for(i in 1:length(out.npn$path)){
			nedges.huge[i] <- sum(out.npn$path[[i]])
		}
		mindist.huge <- which.min(abs(nedges.huge - sum(sim$G)))
		ave.norms <- ave.norms + getNorms(list(out$omega[,,mindist], 
										  as.matrix(out.npn$icov[[mindist.huge]])), Prec)[[2]] 
		ave.f1 <- ave.f1 + getF1(list(out$path[,,mindist],
								as.matrix(out.npn$path[[mindist.huge]])), 
								sim$G)  
		
		p <- dim(X)[2]
		path.huge <- array(NA, dim = c(p, p, length(out.npn$path)))
		omega.huge <- array(NA, dim = c(p, p, length(out.npn$path)))
		for(i in 1:length(out.npn$path)){
			path.huge[,,i] <- as.matrix(out.npn$path[[i]])
			omega.huge[,,i] <- as.matrix(out.npn$icov[[i]])
		}
		edges <- matrix(NA, length(v0s), 2)
		edges.huge <- matrix(NA, length(out.npn$path), 2)
		for(i in 1:length(v0s)){
			edges[i, 1] <- sum(out$path[,,i]) / 2
			edges[i, 2] <- length(intersect(which(out$path[,,i] == 1), which(sim$G == 1))) / 2
		}
		for(i in 1:length(out.npn$path)){
			edges.huge[i, 1] <- sum(path.huge[,,i]) / 2
			edges.huge[i, 2] <- length(intersect(which(path.huge[,,i] == 1), which(sim$G == 1))) / 2
		}
	
		edges <- rbind(c(0,0), edges)
		edges.huge <- rbind(c(0,0), edges.huge)
		edges <- edges[order(edges[, 1]), ]
		edges.huge <- edges.huge[order(edges.huge[, 1]), ]
		# plot(edges[, 1], edges[, 2], type = "l")
		# lines(edges.huge[, 1], edges.huge[, 2], type = "l", col = "red")
		# lines(edges.bd[, 1], edges.bd[, 2], type = "l", col = "blue")

		curve <- approx(edges[,1], edges[,2], xout = 1:(p^2/2), ties = min)$y
		curve[is.na(curve)] <- max(curve, na.rm = TRUE)
		ave.curve <- ave.curve + curve / length(seeds)
		curve2 <- approx(edges.huge[,1], edges.huge[,2], xout = 1:(p^2/2), ties = min)$y
		curve2[is.na(curve2)] <- max(curve2, na.rm = TRUE)
		ave.curve2 <- ave.curve2 + curve2 / length(seeds)
		all.f1[index0, ] <- getF1(list(out$path[,,mindist], 
				   as.matrix(out.npn$path[[mindist.huge]])), sim$G)  
		all.norms[index0,,] <- getNorms(list(out$omega[,,mindist], as.matrix(out.npn$icov[[mindist.huge]])), Prec)[[2]] 
	}
	type <- graph
	name <- paste0(type, ", n = ", n, ", p = ", p)
	print(paste(name, "size:", length(seeds)))
	# plot(1:(p^2/2), ave.curve, type = "l", main = name, xlab = "Total Edge", ylab="TP edges")
	# lines(1:(p^2/2), ave.curve2, lty = 2, col = "blue")
	# legend("bottomright",  c("EMGS", "Glasso"), lty = c(1,2), col = c("black", "blue"))
	ave.f1 <- ave.f1 / count
	ave.norms <- ave.norms/count
	# print(ave.f1)
	# print(ave.norms)
	alllist[[index]] <- list(ave.curve = ave.curve, 
							 ave.curve.npn = ave.curve2,
							 ave.f1 = ave.f1,
							 ave.norms = ave.norms, 
							 all.f1 = all.f1,
							 all.norms = all.norms, 
							 n = n, 
							 p = p, 
							 graph = type)
	names(alllist)[index] <- name
}
library(xml2)
library("knitr")
library("kableExtra")
tab <- matrix(NA, length(ns) * length(simtypes), 3 + 4 * 2)
tab <- data.frame(tab)
tab[, c(1, 2)] <- ""
counter <- 1
for(i in 1:length(ns)){
	for(j in 1){
		tab[counter, 1:2] <- c(ns[i], ps[j]) 
			for(k in 1:length(simtypes)){
			which <- intersect(which(config[,1] == ns[i]), 
							   which(config[,2] == ps[j]))
			which <- intersect(which(config[, 3] == simtypes[k]),
							  which)
			tmp <- alllist[[which]]
			tab[counter, 3] <- simtypes[k]
			tab[counter, 4:5] <- tmp$ave.norms[, 1] 
			tab[counter, 6:7] <- tmp$ave.norms[, 2] 
			tab[counter, 8:9] <- tmp$ave.norms[, 3] 
			tab[counter, 10:11] <- tmp$ave.f1
			counter <- counter + 1
		}
	}
}
tab[is.na(tab[,1]), 1] <- ""
tab[is.na(tab[,2]), 2] <- ""
colnames(tab) <- c("n", "p", "graph", rep(c("EMGS", "NPN"), 4))
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the standardized precision matrix for Graphical model estimation. ", align = "c") %>% add_header_above(c(" ", " ", " ", "M-norm"=2, "S-norm"=2, "F-norm"=2, "F1" = 2)) %>% kable_styling(latex_options = "hold_position")

