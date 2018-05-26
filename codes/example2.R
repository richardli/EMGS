######################################################################
# This script reproduces the example with informative prior
###################################################################### 
library(EMGS)
library(tmvtnorm)
library(mvtnorm)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
library(corrplot)
source("../codes/BDgraph_sim_func.R")

######################################################################
# Simulate three-block precision matrix and data
###################################################################### 
set.seed(12)
p1 <- 10
p2 <- 10
p3 <- 20
p <- p1 + p2 + p3
N <- 100
sim1 <- bdgraph.sim(p = p1, graph = "random", n = 1, type = "Gaussian", prob = 0.8)
sim2 <- bdgraph.sim(p = p2, graph = "random", n = 1, type = "Gaussian", prob = 0.5)
sim3 <- bdgraph.sim(p = p3, graph = "cluster", n = 1, type = "Gaussian", prob = 0.5)
Prec1 <- (sim1$K)
Prec2 <- (sim2$K) 
Prec3 <- (sim3$K)
Prec <- as.matrix(bdiag(Prec1, Prec2, Prec3))
G <- as.matrix(bdiag(matrix(sim1$G, p1, p1), matrix(sim2$G, p2, p2), matrix(sim3$G, p3, p3)))
X <- rtmvnorm(n = N, mean = rep(0, p), sigma = solve(Prec))
v0s <- seq(0.01, 0.2, len = 50)
v1 <- 10
a <- 1
b <- 1
lambda <- 1
groups <- c(rep(0, p1), rep(1, p2), rep(2, p3))

# EMGS without informative prior
out <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.0001, 1e4, verbose = TRUE, copula = FALSE, group = NULL)
# EMGS with informative prior
out1 <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.0001, 1e4, verbose = TRUE, copula = FALSE, group = groups)
# Run glasso with the closest number of edges
out.npn = huge(cov(X), method = "glasso", nlambda=100,lambda.min.ratio = 0.01)
edges.huge <- rep(NA, length(out.npn$path))
for(i in 1:length(out.npn$path)){
	edges.huge[i] <- sum(out.npn$path[[i]]) 
}
mindist2 <- which.min(abs(edges.huge - sum(G)))
lambdas <- seq(out.npn$lambda[mindist2-1], out.npn$lambda[mindist2+1], len = 40) 
out.npn = huge(cov(X), method = "glasso", lambda = lambdas)
edges.huge <- rep(NA, length(out.npn$path))
for(i in 1:length(out.npn$path)){
	edges.huge[i] <- sum(out.npn$path[[i]]) 
}
mindist2 <- which.min(abs(edges.huge - sum(G)))

# Find the iteration with closest number of edges for EMGS
for(i in 1:length(v0s)){
	prob <- getp(v0s[i], v1, out$omega[,,i])
	out$path[,,i][prob >= 0.5] <- 1
	out$path[,,i][prob < 0.5] <- 0
	# Update the posterior probability with Tau. 
	# TODO: need to make this default!
	prob <- getp(v0s[i], v1, out1$omega[,,i], out1$tau[,,i])
	out1$path[,,i][prob >= 0.5] <- 1
	out1$path[,,i][prob < 0.5] <- 0
}
nedges <- apply(out$path, 3, sum)
mindist <- which.min(abs(nedges - sum(G)))
nedges1 <- apply(out1$path, 3, sum)
mindist1 <- which.min(abs(nedges1 - sum(G)))

prec0 <- cov2cor(Prec) * G
prec1 <- cov2cor(out$omega[,,mindist]) 
# prec1[(out$path[,,mindist] <= 0.5)] <- 0
prec1[getp(v0s[mindist], v1, out$omega[,,mindist])<= 0.5] <- 0
diag(prec1) <- 1
prec2 <- cov2cor(out1$omega[,,mindist1])
prec2[getp(v0s[mindist1], v1, out1$omega[,,mindist1], out1$tau[,,mindist1])< 0.5] <- 0

diag(prec2) <- 1
prec3 <- cov2cor(as.matrix(out.npn$icov[[mindist2]]))
diag(prec3) <- 1
diag(prec0) <- 1
prec0[prec0 == 0] <- NA
prec1[prec1 == 0] <- NA
prec2[prec2 == 0] <- NA
prec3[prec3 == 0] <- NA


pdf("../figures/structure.pdf", width = 8, height = 8)
par(mfrow = c(2, 2),oma=c(0,0,2,0)+1)
corrplot(prec3, type = "upper", tl.pos="n", title = "Graphical Lasso", mar = c(0,0,1,0), method="color",  bg = "gray70", na.label.col = "gray70")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70")
corrplot(prec1, type = "upper", tl.pos="n", title = "Exchangeable prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color", bg = "gray70", na.label.col = "gray70")
corrplot(prec2, type = "upper", tl.pos="n", title = "Structured prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70")

par(mar = c(1.8, 0, 0.8, 0)+2)
plot(v0s, out1$tau_compact[1, 1, ], type = "l", col = "red", ylim = c(0, 100), xlab = "v0", ylab = "tau", main = "")
lines(v0s, out1$tau_compact[2, 2, ], col = "orange")
lines(v0s, out1$tau_compact[3, 3, ], col = "blue")
lines(v0s, out1$tau_compact[2, 1, ], col = "darkgreen", lty = 2)
lines(v0s, out1$tau_compact[2, 3, ], col = "darkgreen", lty = 3)
lines(v0s, out1$tau_compact[3, 1, ], col = "darkgreen", lty = 4)
legend("topleft", c("block 1:1","block 2:2","block 3:3", "block 1:2","block 1:3","block 2:3"),
col = c("red", "orange", "blue", rep("darkgreen", 3)), lty = c(1, 1, 1, 2, 3, 4))
abline(h = 1 , col = "gray70", lty = "longdash")
abline(v = v0s[mindist1], col = "gray70", lty = "longdash")
dev.off()
 
