# srun --pty --partition=short --time=3:00:00 --mem-per-cpu=2500 /bin/bash
ns <- c(100, 200, 500)
ps <- c(50, 100, 200)
misses <- c("AR1", "AR2", "random", "cluster")
reps <- 1:100
eval.g <- array(0, dim = c(3, 3, 4, 100))
finish.g <- array(0, dim = c(3, 3, 4, 100))
eval.m <- array(0, dim = c(3, 3, 4, 100))
finish.m <- array(0, dim = c(3, 3, 4, 100))

for(i in 1:3){
for(j in 1:3){
for(k in 1:4){
for(l in 1:100){
	n <- ns[i]
	p <- ps[j]
	graph <- misses[k]
	itr <- reps[l]
	if(file.exists(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))){
		eval.g[i,j,k,l] <- eval.g[i,j,k,l] + 1
		load(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
		finish.g[i,j,k,l] <- finish.g[i,j,k,l] + as.numeric(!is.null(fit$cv.out))

	}
	if(file.exists(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))){
		load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))

		eval.m[i,j,k,l] <- eval.m[i,j,k,l] + 1
		finish.m[i,j,k,l] <- finish.m[i,j,k,l] + as.numeric(!is.null(fit$cv.out))
	}
}
cat(".")
}
}
}


dimnames(eval.g)[[1]] <- dimnames(finish.g)[[1]] <- dimnames(eval.m)[[1]] <- dimnames(finish.m)[[1]] <- paste("n =", ns)
dimnames(eval.g)[[2]] <- dimnames(finish.g)[[2]] <- dimnames(eval.m)[[2]] <- dimnames(finish.m)[[2]] <- paste("p =", ps)
dimnames(eval.g)[[3]] <- dimnames(finish.g)[[3]] <- dimnames(eval.m)[[3]] <- dimnames(finish.m)[[3]] <- misses
g0 <- apply(eval.g, c(1,2,3), sum)
g1 <- apply(finish.g, c(1,2,3), sum)
m0 <- apply(eval.m, c(1,2,3), sum)
m1 <- apply(finish.m, c(1,2,3), sum)

g1/g0
m1/m0
m1-g1
m0-g0