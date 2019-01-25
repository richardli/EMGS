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
ps <- 200
misses <- c("random" )
reps <- 1:100
config <- expand.grid(ns, ps, misses, reps)

allindex <- index  

for(index in allindex){
	print(Sys.time())
	tt <- Sys.time()
	source("simulation-core-1.R")
	print(Sys.time())
	print(Sys.time() - tt)
}
