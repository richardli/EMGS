library(EMGS)
library(ROCR)
library(huge)
library(BDgraph)
source("../codes/functions.R")
index <- 1
# index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
config <- NULL
ns <- c(100, 200, 500)
ps <- 25
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:2
config <- expand.grid(ns, ps, misses, reps)

allindex <- which(c(1:dim(config)[1]) == index)

for(index in allindex){
	source("simcore.R")
	source("evalcore.R")
}