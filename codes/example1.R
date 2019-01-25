######################################################################
# This script reproduces the motivating example
######################################################################
# install.packages("../EMGS_1.0.tar.gz", type = "source", repos = NULL)
set.seed(123)
library(EMGS)
library(huge)
library(tmvtnorm)
# sourceCpp("../EMGS/src/_EMGS.cpp")
N <- 100
M <- 10
Prec <- matrix(0, M, M)
diag(Prec) <- 1
for(i in 2:M){
	Prec[i-1, i] <- 0.5
	Prec[i, i - 1] <- 0.5
}
nedge <- (sum(Prec != 0)-M)/2
diag(Prec) <- 1
# Prec <- Prec * 2
Sigma <- solve(Prec)
X <- rtmvnorm(n = N, mean = rep(0, M), sigma = Sigma)
################################################################
v0s <- seq(0.001, 0.5, len = 40)
v1 <- 100
epsilon <- 1e-6
a <- 2
b <- M
lambda <- 1
omega_path <- array(0, dim = c(length(v0s), M, M))
prob_path <- array(0, dim = c(length(v0s), M, M))
fit <- EMGS(X, v0s, v1, lambda, a, b, epsilon)
for(itr in 1:length(v0s)){
	omega_path[itr, , ] <- tmp <- fit$omega[,,itr]
	tmp[!lower.tri(tmp)] <- NA
	thre <- sort(abs(as.numeric(tmp)), decreasing = TRUE)[nedge]
	prob_path[itr, , ] <- fit$omega[,,itr] >= thre
}
fit.cv <- cv.EMGS(X, v0s, v1, lambda, a, b, epsilon, copula=FALSE, copula_itr = 0)
#########################################################################
glasso.fit <-huge(cov(X), method = "glasso", nlambda=40, lambda.min.ratio = 0.01)
omega_path_glasso0 <- omega_path_glasso <- array(NA, dim = c(M, M, 40))
for(i in 1:dim(omega_path_glasso)[3]){
	omega_path_glasso[,,i] <- as.matrix(glasso.fit$icov[[i]])
	omega_path_glasso0[,,i] <- cov2cor(as.matrix(glasso.fit$icov[[i]]))
	diag(omega_path_glasso0[, , i]) <- NA
}
glasso.cv <- cv.glasso(X = X, rholist = glasso.fit$lambda)

#########################################################################
pdf("../figures/illustration-emgs.pdf", width = 9, height = 4)
par(mfrow = c(1, 2))
omega_path0 <- omega_path
for(i in 1:dim(omega_path0)[1]){
	omega_path0[i, , ] <-  cov2cor(omega_path0[i, , ])
	diag(omega_path0[i, , ]) <- NA
}
plot(v0s, rep(0, length(v0s)), ylim = range(c(-.3, omega_path0, omega_path_glasso0), na.rm = TRUE), col = "white", xlab = "v0", main = "EMGS", ylab = "- Partial correlation")
abline(v = fit.cv$v0.min, lty = 4)
for(i in 2:M){
	for(j in 1:(i-1)){
		if(abs(Prec[i, j]) > 1e-5){
			lines(v0s, omega_path0[, i, j], type = "l", col = "darkblue")
		}
		if(abs(Prec[i, j]) <= 1e-5){
			lines(v0s, omega_path0[, i, j], type = "l", col = "darkgreen", lty = 4)
		}
	}
}
abline(h = 0, col = "red", lty = 2)
abline(h = 0.5, col = "red", lty = 2)
legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("Edges", "Non-edges"))

 
# plot(v0s, rep(0, length(v0s)), ylim = range(c(-.3, omega_path0, omega_path_glasso0), na.rm = TRUE), col = "white", xlab = "v0", main = "EMGS: thresholded\nto have 9 edges", ylab = "- Partial correlation")
# abline(v = fit.cv$v0.min,  lty = 4)
# for(i in 2:M){
# 	for(j in 1:(i-1)){
# 		if(abs(Prec[i, j]) > 1e-5){
# 			lines(v0s, omega_path0[, i, j] * (prob_path[, i, j] > 0.5), type = "l", col = "darkblue")
# 		}
# 		if(abs(Prec[i, j]) <= 1e-5){
# 			lines(v0s, omega_path0[, i, j] * (prob_path[, i, j] > 0.5), type = "l", col = "darkgreen", lty = 4)
# 		}
# 	}
# }
# abline(h = 0, col = "red", lty = 2)
# abline(h = 0.5, col = "red", lty = 2)
# legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("Edges", "Non-edges"))

plot(glasso.fit$lambda, rep(0, length(glasso.fit$lambda)), ylim = range(c(-.3, omega_path0, omega_path_glasso0), na.rm = TRUE), col = "white", xlab = "rho", main = "Graphical Lasso", ylab = "- Partial correlation")
abline(v = glasso.cv$rho.min, lty = 4)
# abline(v = glasso.cv$rho.1se, lty = 4)
for(i in 2:M){
	for(j in 1:(i-1)){
		if(Prec[i, j] != 0){
		lines(glasso.fit$lambda, omega_path_glasso0[i, j, ], type = "l", col = "darkblue")			
		}
		if(Prec[i, j] == 0){
		lines(glasso.fit$lambda, omega_path_glasso0[i, j, ], type = "l", col = "darkgreen", lty = 4)			
		}
	}
}
abline(h = 0, col = "red", lty = 2)
abline(h = 0.5, col = "red", lty = 2)
legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("Edges", "Non-edges"))
dev.off()