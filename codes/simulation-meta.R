# srun --pty --partition=short --time=3:00:00 --mem-per-cpu=2500 /bin/bash
# srun --pty --partition=build --time=1:00:00 --mem-per-cpu=2500 /bin/bash
# module load R
# R
# install.packages("../EMGS_0.0.1.tar.gz", repos = NULL, type = "source")

#!/bin/bash
#SBATCH --job-name sim2 # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition short     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-05:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=3000     # Memory limit for each tasks (in MB)
#SBATCH -o out/Rt_%j.out    # File to which STDOUT will be written
#SBATCH -e out/Rt_%j.err    # File to which STDERR will be written

# module load R
# Rscript runbatch2.R $SLURM_ARRAY_TASK_ID > ../experiments/log/g-$SLURM_ARRAY_TASK_ID

# sbatch --array=1-500 runsim2.sbatch
# sbatch --array=1-500 runsim.sbatch
# sbatch --array=500-1000 runsim2.sbatch
# sbatch --array=500-1000 runsim.sbatch

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
index <- 1
index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
config <- NULL
ns <- c(100, 200, 500)
ps <- c(25, 50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:100
config <- expand.grid(ns, ps, misses, reps)

allindex <- which(c(1:dim(config)[1]) %% 1000 == index)

for(index in allindex){
	source("simulation-core.R")
}


# for(index in allindex){
# args_in <- as.character(config[index, ])

# p = as.numeric(args_in[2])
# n = as.numeric(args_in[1])
# graph = args_in[3]
# itr = as.numeric(args_in[4])
# set.seed(itr * 123456)
# load( paste("../data/simG/graphsim", n, p, graph, itr, ".rda", sep="_"))
# X <- sim$data
# v0s <- seq(0.001, 0.2, len = 40)
# v1 <- 1000
# a <- 1
# b <- 1
# lambda <- 1
# Prec <- sim$K  

# fit.cv <- cv.EMGS(X, v0s, v1, lambda, a, b, copula = FALSE, copula_itr = 0)
# fit.glasso.cv <- cv.glasso(X)
# fit <- list(out = fit.cv, out.npn = fit.glasso.cv, bdgraph.obj = NULL)
# save(fit, file = paste("../data/simGresults/graphfitcv", n, p, graph, itr, ".rda", sep="_"))
# }


# ############
# bdgraph.obj <- bdgraph(data = sim, iter = 1e4, method = "gcgm")
# fit <- list(out = out, out.npn = out.npn, bdgraph.obj = bdgraph.obj)
# save(fit, file = paste("../data/sim/graphfit", n, p, graph, itr, ".rda", sep="_"))



################################################################
#  Select by sparsity
################################################################
remove(list = ls())
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
config <- NULL
ns <- c(50, 100, 500)
ps <- c(10, 50)
misses <- c("AR1", "AR2", "random", "cluster")
config <- expand.grid(ns, ps, misses)
v0s <- seq(0.001, 0.2, len = 40)

files <- list.files("../data/simGresults")
alllist <- vector("list", dim(config)[1])
for(index in 1:dim(config)[1]){
	args_in <- as.character(config[index, ])

	p = as.numeric(args_in[2])
	n = as.numeric(args_in[1])
	graph = args_in[3]
	itr = as.numeric(args_in[4])
	tmp <- gsub(paste("graphfit", n, p, graph, "", sep="_"), "!", files)
	tmp <- gsub("_.rda", "!", tmp)
	tmp <- unlist(strsplit(tmp, "!"))
	# force the non-rep number to be NA
	seeds <- sort(unique(as.numeric(tmp)))
	ave.curve <- ave.curve2 <- rep(0, p^2/2)
	ave.norms <- matrix(0, 2, 4)
	ave.f1 <- rep(0, 2)
	all.norms <- array(NA, c(length(seeds), 2, 4))
	all.f1 <- matrix(NA, length(seeds), 2)
	count <- 0
	for(index0 in 1:length(seeds)){
		itr <- seeds[index0]
		load( paste("../data/simGresults/graphfit", n, p, graph, itr, ".rda", sep="_"))
		load( paste("../data/simG/graphsim", n, p, graph, itr, ".rda", sep="_"))
		if(sum(sim$G) < 1){next}
		count <- count + 1

		Prec <- sim$K  
		X <- sim$data
		out <- fit$out
		out.npn <- fit$out.npn
		nedges <- apply(out$path, 3, sum)
		mindist <- which.min(abs(nedges - sum(sim$G)))
		getNorms(list(out$omega[,,mindist]), Prec)[[1]]
		getF1(list(out$path[,,mindist]), sim$G)
		# path.plot(v0s, out, sim$G, thres = 0.5, normalize=TRUE, main = paste("N =",  dim(X)[1], "M =", dim(X)[2]))

		nedges.huge <- rep(0, length(out.npn$path))
		for(i in 1:length(out.npn$path)){
			nedges.huge[i] <- sum(out.npn$path[[i]])
		}
		mindist.huge <- which.min(abs(nedges.huge - sum(sim$G)))
		ave.norms <- ave.norms + getNorms(list(out$omega[,,mindist], as.matrix(out.npn$icov[[mindist.huge]])), Prec)[[2]] 
		ave.f1 <- ave.f1 + getF1(list(out$path[,,mindist], 
				   as.matrix(out.npn$path[[mindist.huge]])), sim$G)  
		
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
		# cutoff <- sort(unique(as.numeric(bdgraph.obj$p_links)))
		# edges.bd <- matrix(NA, length(cutoff), 2)
		# for(i in 1:length(cutoff)){
		# 	# upper diagonal here
		# 	edges.bd[i, 1] <- sum(bdgraph.obj$p_links > cutoff[i]) 
		# 	edges.bd[i, 2] <- length(intersect(which(bdgraph.obj$p_links > cutoff[i]), which(sim$G == 1))) 
		# }
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
	type <-  misses[as.numeric(graph)]
	name <- paste0(type, ", n = ", n, ", p = ", p)
	print(paste(name, "size:", length(seeds)))
	plot(1:(p^2/2), ave.curve, type = "l", main = name)
	lines(1:(p^2/2), ave.curve2, lty = 2, col = "blue")
	ave.f1 <- ave.f1 / count
	ave.norms <- ave.norms/count
	print(ave.f1)
	print(ave.norms)
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

save(alllist, file = "../data/simG0.rda")


load("../data/simG0.rda")
# install.packages("xml2")
library(xml2)
library("knitr")
library("kableExtra")
config <- NULL
ns <- c(50, 100, 500)
ps <- c(10, 50)
misses <- c("AR1", "AR2", "random", "cluster")
config <- expand.grid(ns, ps, misses)

tab <- matrix(NA, length(ns) * length(misses), 3 + 4 * 2)
tab <- data.frame(tab)
tab[, c(1, 2)] <- ""
counter <- 1
for(i in 1:length(ns)){
	for(j in 1){
		tab[counter, 1:2] <- c(ns[i], ps[j]) 
			for(k in 1:length(misses)){
			which <- intersect(which(config[,1] == ns[i]), 
							   which(config[,2] == ps[j]))
			which <- intersect(which(config[, 3] == misses[k]),
							  which)
			tmp <- alllist[[which]]
			tab[counter, 3] <- misses[k]
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


colnames(tab) <- c("n", "p", "graph", rep(c("EMGS", "Glasso"), 4))
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the standardized precision matrix for Graphical model estimation with small $p$. The final graphs are chosen so that the sparsity level is closest to the truth. ", align = "c") %>% add_header_above(c(" ", " ", " ", "M-norm"=2, "S-norm"=2, "F-norm"=2, "F1" = 2)) %>% kable_styling(latex_options = "hold_position")



tab <- matrix(NA, length(ns) * length(misses), 3 + 4 * 2)
tab <- data.frame(tab)
tab[, c(1, 2)] <- ""
counter <- 1
for(i in 1:length(ns)){
	for(j in 2){
		tab[counter, 1:2] <- c(ns[i], ps[j]) 
			for(k in 1:length(misses)){
			which <- intersect(which(config[,1] == ns[i]), 
							   which(config[,2] == ps[j]))
			which <- intersect(which(config[, 3] == misses[k]),
							  which)
			tmp <- alllist[[which]]
			tab[counter, 3] <- misses[k]
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


colnames(tab) <- c("n", "p", "graph", rep(c("EMGS", "Glasso"), 4))
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the standardized precision matrix for Graphical model estimation with moderate $p$. The final graphs are chosen so that the sparsity level is closest to the truth.", align = "c") %>% add_header_above(c(" ", " ", " ", "M-norm"=2, "S-norm"=2, "F-norm"=2, "F1" = 2)) %>% kable_styling(latex_options = "hold_position")



toplot <- data.frame(mean=NULL, lower=NULL, upper=NULL, method=NULL, name=NULL,n=NULL, p=NULL,type=NULL,name=NULL)
for(i in 1:length(alllist)){

	for(j in 1:4){
		mean <- apply(alllist[[i]]$all.norms[,,j], 2, mean,na.rm=T)
		lower <- apply(alllist[[i]]$all.norms[,,j], 2, quantile, 0.025,na.rm=T)
		upper <- apply(alllist[[i]]$all.norms[,,j], 2, quantile, 0.975,na.rm=T)

		tmp <- data.frame(mean=mean, lower = lower, upper = upper, type = alllist[[i]]$graph, n = alllist[[i]]$n,p = alllist[[i]]$p, method = c("EMGS", "GLasso"), name = colnames(alllist[[i]]$ave.norms)[j])
		toplot <- rbind(toplot, tmp)
	}
	mean <- apply(alllist[[i]]$all.f1, 2, mean,na.rm=T)
	lower <- apply(alllist[[i]]$all.f1, 2, quantile, 0.025,na.rm=T)
	upper <- apply(alllist[[i]]$all.f1, 2, quantile, 0.975,na.rm=T)

	tmp <- data.frame(mean=mean, lower = lower, upper = upper, type = alllist[[i]]$graph, n = alllist[[i]]$n,p = alllist[[i]]$p, method = c("EMGS", "GLasso"), name = "F-1")
	toplot <- rbind(toplot, tmp)
}

library(ggplot2)

pdf("../figures/g50.pdf", width = 8, height = 5)
toplot$n <- factor(toplot$n)
g <- ggplot(data = subset(toplot, p == 50 & name %in% c("M norm", "F norm", "F-1")), aes(x=n,y=mean,group=method,color=method))
g <- g + geom_point(position=position_dodge(width=0.3), size=1)
g <- g + geom_line(position=position_dodge(width=0.3), size=0.8)
g <- g + geom_errorbar(aes(ymin=lower, ymax=upper, color = method), width=.1, position=position_dodge(width=0.3), size = 0.5) 
g <- g + facet_grid(name~type, scale="free")
g <- g + ggtitle("Gaussian graphical model: p = 50") + ylab("") + xlab("Sample size")
g
dev.off()




################################################################
#  Select by CV
################################################################
remove(list = ls())
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
config <- NULL
ns <- c(50, 100, 500)
ps <- c(10, 50)
misses <- c("AR1", "AR2", "random", "cluster")
config <- expand.grid(ns, ps, misses)
v0s <- seq(0.001, 0.2, len = 40)

files <- list.files("../data/simGresults")
alllist <- vector("list", dim(config)[1])
for(index in 1:dim(config)[1]){
	args_in <- as.character(config[index, ])

	p = as.numeric(args_in[2])
	n = as.numeric(args_in[1])
	graph = args_in[3]
	itr = as.numeric(args_in[4])
	tmp <- gsub(paste("graphfitcv", n, p, graph, "", sep="_"), "!", files)
	tmp <- gsub("_.rda", "!", tmp)
	tmp <- unlist(strsplit(tmp, "!"))
	# force the non-rep number to be NA
	seeds <- sort(unique(as.numeric(tmp)))
	if(length(seeds) == 0){
		next
	}
	ave.curve <- ave.curve2 <- rep(0, p^2/2)
	ave.norms <- matrix(0, 3, 4)
	ave.f1 <- rep(0, 3)

	count <- 0
	for(itr in seeds){
		load( paste("../data/simGresults/graphfitcv", n, p, graph, itr, ".rda", sep="_"))
		load( paste("../data/simG/graphsim", n, p, graph, itr, ".rda", sep="_"))
		if(sum(sim$G) < 1){next}
		count <- count + 1

		Prec <- sim$K  
		X <- sim$data
		out <- fit$out
		out.npn <- fit$out.npn

		if(is.null(fit$out.bd)){
			bdgraph.obj <- bdgraph(data = sim, iter = 1e4)
			fit$out.bd <- bdgraph.obj
			save(fit, file = paste("../data/simGresults/graphfitcv", n, p, graph, itr, ".rda", sep="_"))
		}

		ave.norms <- ave.norms + getNorms(list(
			out$fit.1se$omega[,,1], 
			as.matrix(out.npn$fit$wi[,,1]),
			fit$out.bd$K_hat
			), Prec)[[2]] 

		path3 <- fit$out.bd$p_links 
		path3 <- path3 + t(path3)
		path3 <- path3 > 0.5
		ave.f1 <- ave.f1 + getF1(list(
			out$fit.1se$path[,,1], 
			as.matrix(out.npn$fit$wi[,,1]),
			path3
			), sim$G)  
	}
	type <-  misses[as.numeric(graph)]
	name <- paste0(type, ", n = ", n, ", p = ", p)
	print(name)
	ave.f1 <- ave.f1 / count
	ave.norms <- ave.norms/count
	print(ave.f1)
	print(ave.norms)
	alllist[[index]] <- list(ave.f1 = ave.f1,
							 ave.norms = ave.norms,
							 n = n, 
							 p = p, 
							 graph = type)
	names(alllist)[index] <- name
}




library(xml2)
library("knitr")
library("kableExtra")
config <- NULL
ns <- c(50, 100, 500)
ps <- c(10, 50)
misses <- c("AR1", "AR2", "random", "cluster")
config <- expand.grid(ns, ps, misses)

tab <- matrix(NA, length(n) * length(misses), 3 + 4 * 3)
tab <- data.frame(tab)
tab[, c(1, 2)] <- ""
counter <- 1
for(i in 1:length(ns)){
	for(j in 1){
		tab[counter, 1:2] <- c(ns[i], ps[j]) 
			for(k in 1:length(misses)){
			which <- intersect(which(config[,1] == ns[i]), 
							   which(config[,2] == ps[j]))
			which <- intersect(which(config[, 3] == misses[k]),
							  which)
			tmp <- alllist[[which]]
			tab[counter, 3] <- misses[k]
			if(!is.null(alllist[[which]])){
				tab[counter, 4:6] <- tmp$ave.norms[, 1] 
				tab[counter, 7:9] <- tmp$ave.norms[, 2] 
				tab[counter, 10:12] <- tmp$ave.norms[, 3] 
				tab[counter, 13:15] <- tmp$ave.f1				
			}
			counter <- counter + 1
		}
	}
}
tab[is.na(tab[,1]), 1] <- ""
tab[is.na(tab[,2]), 2] <- ""


colnames(tab) <- c("n", "p", "graph", rep(c("EMGS", "Glasso", "BDMCMC"), 4))
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the standardized precision matrix for Graphical model estimation. Graph selected by cross-validation", align = "c") %>% add_header_above(c(" ", " ", " ", "M-norm"=2, "S-norm"=2, "F-norm"=2, "F1" = 2)) %>% kable_styling(latex_options = "hold_position")


tab <- matrix(NA, length(n) * length(misses), 3 + 4 * 3)
tab <- data.frame(tab)
tab[, c(1, 2)] <- ""
counter <- 1
for(i in 1:length(ns)){
	for(j in 2){
		tab[counter, 1:2] <- c(ns[i], ps[j]) 
			for(k in 1:length(misses)){
			which <- intersect(which(config[,1] == ns[i]), 
							   which(config[,2] == ps[j]))
			which <- intersect(which(config[, 3] == misses[k]),
							  which)
			tmp <- alllist[[which]]
			tab[counter, 3] <- misses[k]
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


colnames(tab) <- c("n", "p", "graph", rep(c("EMGS", "Glasso", "BDMCMC"), 4))
kable(tab, format = "latex", digits = 2, booktabs = TRUE, caption = "Comparing estimation of the standardized precision matrix for Graphical model estimation with moderate $p$. The final graphs are chosen by cross-validation.", align = "c") %>% add_header_above(c(" ", " ", " ", "M-norm"=2, "S-norm"=2, "F-norm"=2, "F1" = 2)) %>% kable_styling(latex_options = "hold_position")

