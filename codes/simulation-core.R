args_in <- config[index, ]
p = as.numeric(args_in[2])
n = as.numeric(args_in[1])
graph = as.character(unlist(args_in[3]))
itr = as.numeric(args_in[4])
print(paste("Check", n, p, graph, itr, sep="_"))

###################################################################
## Gaussian case
run <- FALSE
if(file.exists(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))){
	load(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
	if(is.null(fit$cv.out)){
		run <- TRUE
		print(paste("Refit", n, p, graph, itr, sep="_"))
	}
}else{
	run <- TRUE
	print(paste("Refit", n, p, graph, itr, sep="_"))
}
if(run){
	set.seed(itr * 123456)
	load( paste("../data/simG/graphsim", n, p, graph, itr, "gaussian.rda", sep="_"))
	X <- sim$data
	v0s <- seq(0.01, 1, len = 40)
	v1 <- 100
	a <- 1
	b <- 1
	lambda <- 1
	Prec <- sim$K  
	# Glasso
	out.npn = huge(X, method = "glasso", nlambda=40)
	cv.out.npn = cv.glasso(X, rholist = out.npn$lambda, K=5)

	fit <- list(out = NULL, out.npn = out.npn, 
				cv.out = NULL, cv.out.npn = cv.out.npn,
				Prec = Prec, G = sim$G)
	save(fit, file = paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))

	# EMGS
	out <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.001, 1e4, verbose = TRUE, copula = FALSE, copula_itr = 30)
	cv.out <- cv.EMGS(as.matrix(X), v0s, v1, lambda, a, b, epsilon=0.001, K=5, maxitr=1e4,  copula = FALSE, copula_itr=30)
	fit <- list(out = out, out.npn = out.npn, 
				cv.out = cv.out, cv.out.npn = cv.out.npn,
				Prec = Prec, G = sim$G)
	save(fit, file = paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
}




###################################################################
## Mixed case
run <- FALSE
if(file.exists(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))){
	load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
	if(is.null(fit$cv.out)){
		run <- TRUE
		print(paste("Mixed Refit", n, p, graph, itr, sep="_"))
	}
}else{
	run <- TRUE
	print(paste("Mixed Refit", n, p, graph, itr, sep="_"))
}
if(run){
	set.seed(itr * 123456)
	load( paste("../data/simM/graphsim", n, p, graph, itr, "mixed.rda", sep="_"))
	X <- sim$data
	v0s <- seq(0.01, 1, len = 40)
	v1 <- 100
	a <- 1
	b <- 1
	lambda <- 1
	Prec <- sim$K  

	## Glasso
	transform <- huge.npn(X)
	out.npn = huge(transform, method = "glasso", nlambda=40)
	cv.out.npn = cv.glasso(transform, rholist = out.npn$lambda, K=5)
	 
	fit <- list(out = NULL, out.npn = out.npn, 
				cv.out = NULL, cv.out.npn = cv.out.npn,  
				Prec = Prec, G = sim$G)
	save(fit, file = paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))

	## EMGS
	out <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.01, 1e4, verbose = TRUE, copula = TRUE, copula_itr = 30)
	cv.out <- cv.EMGS(as.matrix(X), v0s, v1, lambda, a, b, epsilon=0.01, K=5, maxitr=1e4,  copula = TRUE, copula_itr=30)
	fit <- list(out = out, out.npn = out.npn, 
				cv.out = cv.out, cv.out.npn = cv.out.npn,
				Prec = Prec, G = sim$G)
	save(fit, file = paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
}
