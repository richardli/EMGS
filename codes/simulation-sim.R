library(BDgraph)
source("BDgraph_sim_func.R")
ns <- c(100, 200, 500)
ps <- c(25, 50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
config <- expand.grid(ns, ps, misses)
Nitr <- 1

set.seed(123456)
for(itr in 1:Nitr){
	cat(".")
	for(index in 1:dim(config)[1]){
		args_in <- config[index, ]
		p = as.numeric(args_in[2])
		n = as.numeric(args_in[1])
		graph = as.character(unlist(args_in[3]))
		##########################################################
		sim <- bdgraph.sim(p = p, graph = graph, n = n, type = "Gaussian", prob = 0.2)
		sim$sigma <- cov2cor(solve(sim$K))
		sim$K <- solve(sim$sigma)
		sim$data <- mvtnorm::rmvnorm(n=n, mean = rep(0, p), sigma = sim$sigma)
		save(sim, file = paste("../data/simG/graphsim", n, p, graph, itr, "gaussian.rda", sep="_"))
		##################################################
        d <- sim$data
        ps = floor(p/2)
        col_number <- c(1:ps)
        prob <- stats::pnorm(d[, col_number])
        d[, col_number] <- stats::qpois(p = prob, lambda = 10)
        col_number <- c((ps + 1):p)
        prob <- stats::pnorm(d[, col_number])
        d[, col_number] <- stats::qpois(p = prob, lambda = 2)
        # col_number <- c((2 * ps + 1):p)
        # prob <- stats::pnorm(d[, col_number])
        # d[, col_number] <- stats::qexp(p = prob, rate = 10)
        sim$data <- d
        for (j in 1:p) {
        	 sim$data[, j] <- sim$data[, j] - mean(sim$data[, j])
        }
		save(sim, file = paste("../data/simM/graphsim", n, p, graph, itr, "mixed.rda", sep="_"))
	}
}