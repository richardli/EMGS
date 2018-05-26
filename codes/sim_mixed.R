##########################################################################
## Install package
##########################################################################
# install.packages("../EMGS_1.0.tar.gz", repos = NULL, type = "source")

library(EMGS)
library(tmvtnorm)
library(mvtnorm)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
source("../codes/BDgraph_sim_func.R")
##########################################################################
## Simulate and save data files
##########################################################################
config <- NULL
ns <- c(50)
ps <- c(50)
simtypes <- c("AR1", "AR2", "random", "cluster")
Nsim <- 10
config <- expand.grid(ns, ps, simtypes, stringsAsFactors = FALSE)

set.seed(123456)
for(itr in 1:Nsim){
	for(index in 1:dim(config)[1]){
		args_in <- config[index, ]
		p = as.numeric(args_in[2])
		n = as.numeric(args_in[1])
		graph = as.character(args_in[3])
		sim <- bdgraph.sim(p = p, graph = graph, n = n, type = "mixed", prob = 0.2)
		save(sim, file = paste("../data/sim/graphsim", n, p, graph, itr, "mixed.rda", sep="_"))
	}
}

##########################################################################
## Fit models for one configuration over 10 replications
##########################################################################
config <- NULL
ns <- c(50)
ps <- c(50)
simtypes <- "random"  # c("AR1", "AR2", "random", "cluster")
reps <- 1:10
config <- expand.grid(ns, ps, simtypes, reps, stringsAsFactors = FALSE)
allindex <- 1:dim(config)[1]

for(index in allindex){
args_in <- as.character(config[index, ])

p = as.numeric(args_in[2])
n = as.numeric(args_in[1])
graph = args_in[3]
itr = as.numeric(args_in[4])
set.seed(itr * 123456)
load( paste("../data/sim/graphsim", n, p, graph, itr, "mixed.rda", sep="_"))
X <- sim$data
v0s <- seq(0.001, 0.2, len = 40)
v1 <- 1000
a <- 1
b <- 1
lambda <- 1
Prec <- sim$K  
out <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.01, 1e4, verbose = TRUE, copula = TRUE, copula_itr = p)
transform <- huge.npn(X)
out.npn = huge(transform, method = "glasso", nlambda=40,lambda.min.ratio = 0.4)
# bdgraph.obj <- bdgraph(data = sim, iter = 1e4, method = "gcgm")
# fit <- list(out = out, out.npn = out.npn, bdgraph.obj = bdgraph.obj)
save(fit, file = paste("../data/fit/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
}


