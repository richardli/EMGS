library(EMGS)
library(mvtnorm)
library(tmvtnorm)
library(ROCR)
library(huge)
library(BDgraph)
library(truncnorm)
library(corrplot)
source("../codes/functions.R")

######################################################################
# Simulate three-block precision matrix and data
###################################################################### 
set.seed(123)
p1 <- 20
p2 <- 20
p3 <- 20
p <- p1 + p2 + p3
N <- 200
sim1 <- bdgraph.sim(p = p1, graph = "random", n = 1, type = "Gaussian", prob = 0.4)
sim2 <- bdgraph.sim(p = p2, graph = "random", n = 1, type = "Gaussian", prob = 0.4)
sim3 <- bdgraph.sim(p = p3, graph = "random", n = 1, type = "Gaussian", prob = 0.4)
Prec1 <- cov2cor(sim1$K) 
Prec2 <- cov2cor(sim2$K) * .7  
diag(Prec2) <- 1
# diag(Prec2) <- diag(sim2$K)
Prec3 <- cov2cor(sim3$K) * .4
diag(Prec3) <- 1
# diag(Prec3) <- diag(sim3$K)
Prec <- as.matrix(bdiag(Prec1, Prec2, Prec3))
Cov <- solve(Prec)
Prec <- solve(cov2cor(Cov))
# corrplot(Prec, is.corr=F, method = "color")
G <- as.matrix(bdiag(matrix(sim1$G, p1, p1), matrix(sim2$G, p2, p2), matrix(sim3$G, p3, p3)))
X <- rtmvnorm(n = N, mean = rep(0, p), sigma = solve(Prec))
v0s <- seq(0.01, 1, len = 40)
v1 <- 100
a <- 1
b <- 10
lambda <- 1
groups <- c(rep(0, p1), rep(1, p2), rep(2, p3))
nedge <- sum(G)/2

# EMGS without informative prior
out <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.0001, 1e4, verbose = TRUE, copula = FALSE, group = NULL)
out.cv <- cv.EMGS(as.matrix(X), v0s, v1, lambda, a, b, epsilon=0.0001, K=5, maxitr=1e4,  copula = FALSE, group = NULL, copula_itr=0)

# EMGS with informative prior
out1 <- EMGS(as.matrix(X), v0s, v1, lambda, a, b, 0.0001, 1e4, verbose = TRUE, copula = FALSE, group = groups)
out1.cv <- cv.EMGS(as.matrix(X), v0s, v1, lambda, a, b, epsilon=0.0001, K=5, maxitr=1e4,  copula = FALSE, group = groups, copula_itr=0, a_tau=1, b_tau = 1)

# Run glasso
out.npn = huge(cov(X), method = "glasso", nlambda=100,lambda.min.ratio = 0.01)
out.npn.cv = cv.glasso(X, rholist = out.npn$lambda, K=5)

if(FALSE){
	print(c(out.cv$v0.min, out.cv$v0.1se))
	print(c(out1.cv$v0.min, out1.cv$v0.1se))
	print(c(out.npn.cv$rho.min, out.npn.cv$rho.1se))
	plot(v0s, apply(out1.cv$ll, 1, mean), type = "l")
	lines(v0s, apply(out.cv$ll, 1, mean), col="gray")
}

prec0 <- -cov2cor(Prec) * G
diag(prec0) <- 1

prec1 <- -cov2cor(out.cv$fit.min$omega[,,1])
tmp <- out.cv$fit.min$omega[,,1] 
tmp[!lower.tri(tmp)] <- NA
thre <- sort(abs(as.numeric(tmp)), decreasing = TRUE)[nedge]
prec1[abs(out.cv$fit.min$omega[,,1]) < thre] <- NA
diag(prec1) <- 1

# getThres <- function(p, v1, v0, tau){
# 	return(sqrt(
# 			(log(p/(1-p)) + log(v1/v0)) * 2*v0^2*v1^2/(v1^2-v0^2) / tau
# 		)
# 	)
# }
# thres1 <- getThres(0.5, out1.cv$fit.min$v0, out1.cv$fit.min$v1, out1.cv$fit.min$tau[,,1])
# thres <- getThres(0.9, out.cv$fit.min$v0, out.cv$fit.min$v1, matrix(1, p, p))
# thres[1]
# corrplot(out.cv$fit.min$omega[,,1], is.corr=F, method="color")
# corrplot(thres, is.corr=F, method="color")
# tmp <- out.cv$fit.min$omega[,,1]
# tmp[abs(tmp) < thres] <- NA
# corrplot(as.matrix(tmp), is.corr=F, method="color")


prec2 <- -cov2cor(out1.cv$fit.min$omega[,,1])
tmp <- out1.cv$fit.min$omega[,,1] 
tmp[!lower.tri(tmp)] <- NA
thre <- sort(abs(as.numeric(tmp)), decreasing = TRUE)[nedge]
prec2[abs(out1.cv$fit.min$omega[,,1]) < thre] <- NA
diag(prec2) <- 1


prec3 <- -cov2cor(as.matrix(out.npn.cv$fit.min$icov[[1]]))
tmp <- as.matrix(out.npn.cv$fit.min$icov[[1]]) 
tmp[!lower.tri(tmp)] <- NA
thre <- sort(abs(as.numeric(tmp)), decreasing = TRUE)[nedge]
prec3[abs(as.matrix(out.npn.cv$fit.min$icov[[1]])) < thre] <- NA
diag(prec3) <- 1

prec0[prec0 == 0] <- NA
prec1[prec1 == 0] <- NA
prec2[prec2 == 0] <- NA
prec3[prec3 == 0] <- NA


pdf("../figures/structure.pdf", width = 8, height = 8)
par(mfrow = c(2, 2),oma=c(0,0,2,0))
corrplot(prec3, type = "upper", tl.pos="n", title = "Graphical Lasso", mar = c(0,0,1,0), method="color",  bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec1, type = "upper", tl.pos="n", title = "Exchangeable prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color", bg = "gray70", na.label.col = "gray70", na.label = " ")
corrplot(prec2, type = "upper", tl.pos="n", title = "Structured prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70", na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70", na.label = " ")

par(mar = c(1.8, 0, 0.8, 0)+2)
if(FALSE){
plot(v0s, out1$tau_compact[1, 1, ], type = "l", col = "red", ylim = range(c(0, 20)), xlab = "v0", ylab = "tau", main = "")
lines(v0s, out1$tau_compact[2, 2, ], col = "orange")
lines(v0s, out1$tau_compact[3, 3, ], col = "blue")
lines(v0s, out1$tau_compact[2, 1, ], col = "darkgreen", lty = 2)
lines(v0s, out1$tau_compact[2, 3, ], col = "darkgreen", lty = 3)
lines(v0s, out1$tau_compact[3, 1, ], col = "darkgreen", lty = 4)
}

plot(v0s, 1/out1$tau_compact[1, 1, ]^.5, type = "l", col = "red", ylim = range(1/out1$tau_compact), xlab = "v0", ylab = "tau", main = "")
lines(v0s, 1/out1$tau_compact[2, 2, ]^.5, col = "orange")
lines(v0s, 1/out1$tau_compact[3, 3, ]^.5, col = "blue")
lines(v0s, 1/out1$tau_compact[2, 1, ]^.5, col = "darkgreen", lty = 2)
lines(v0s, 1/out1$tau_compact[2, 3, ]^.5, col = "darkgreen", lty = 3)
lines(v0s, 1/out1$tau_compact[3, 1, ]^.5, col = "darkgreen", lty = 4)

legend("topright", c("block 1:1", "block 2:2", "block 3:3", "block 1:2", "block 1:3", "block 2:3"), col = c("red", "orange", "blue", rep("darkgreen", 3)), lty = c(1, 1, 1, 2, 3, 4))
abline(v = out1.cv$v0.min, col = "gray70", lty = "longdash")
dev.off()
 


pdf("../figures/structure2.pdf", width = 12, height = 3)
par(mfrow = c(1, 4),oma=c(0,0,2,0))
corrplot(prec3, type = "upper", tl.pos="n", title = "Graphical Lasso", mar = c(0,0,1,0), method="color",  bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec1, type = "upper", tl.pos="n", title = "Exchangeable prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70",na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color", bg = "gray70", na.label.col = "gray70", na.label = " ")
corrplot(prec2, type = "upper", tl.pos="n", title = "Structured prior", mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70", na.label = " ")
corrplot(prec0, type = "lower", add=TRUE, tl.pos="n", cl.pos="n", method="color",  bg = "gray70", na.label.col = "gray70", na.label = " ")

par(mar = c(1.8, 0, 0.8, 0)+2)
if(FALSE){
plot(v0s, out1$tau_compact[1, 1, ], type = "l", col = "red", ylim = range(c(0, 20)), xlab = "v0", ylab = "tau", main = "")
lines(v0s, out1$tau_compact[2, 2, ], col = "orange")
lines(v0s, out1$tau_compact[3, 3, ], col = "blue")
lines(v0s, out1$tau_compact[2, 1, ], col = "darkgreen", lty = 2)
lines(v0s, out1$tau_compact[2, 3, ], col = "darkgreen", lty = 3)
lines(v0s, out1$tau_compact[3, 1, ], col = "darkgreen", lty = 4)
}

plot(v0s, 1/out1$tau_compact[1, 1, ]^.5, type = "l", col = "red", ylim = range(1/out1$tau_compact), xlab = "v0", ylab = "tau", main = "scaling factor")
lines(v0s, 1/out1$tau_compact[2, 2, ]^.5, col = "orange")
lines(v0s, 1/out1$tau_compact[3, 3, ]^.5, col = "blue")
lines(v0s, 1/out1$tau_compact[2, 1, ]^.5, col = "darkgreen", lty = 2)
lines(v0s, 1/out1$tau_compact[2, 3, ]^.5, col = "darkgreen", lty = 3)
lines(v0s, 1/out1$tau_compact[3, 1, ]^.5, col = "darkgreen", lty = 4)

legend("topright", c("block 1:1", "block 2:2", "block 3:3", "block 1:2", "block 1:3", "block 2:3"), col = c("red", "orange", "blue", rep("darkgreen", 3)), lty = c(1, 1, 1, 2, 3, 4))
abline(v = out1.cv$v0.min, col = "gray70", lty = "longdash")
dev.off()
 







