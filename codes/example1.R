######################################################################
# This script reproduces the motivating example
######################################################################
library(EMGS)
library(huge)
library(tmvtnorm)

N <- 100
M <- 10
Prec <- matrix(0, M, M)
diag(Prec) <- 1
for(i in 2:M){
	Prec[i-1, i] <- 0.5
	Prec[i, i - 1] <- 0.5
}
diag(Prec) <- 1
Sigma <- solve(Prec)
X <- rtmvnorm(n = N, mean = rep(0, M), sigma = Sigma)
################################################################
v0s <- seq(0.001, 0.15, len = 40)
v1 <- 10000
epsilon <- 1e-6
a <- 1
b <- 1
lambda <- 10
omega_path <- array(0, dim = c(length(v0s), M, M))
prob_path <- array(0, dim = c(length(v0s), M, M))
for(itr in 1:length(v0s)){
	################################ 
	v0 <- v0s[itr]
	fit <- EMGS(X, v0, v1, lambda, a, b, epsilon)
	omega_path[itr, , ] <- fit$omega
	prob_path[itr, , ] <- fit$EZ
}
#########################################################################
glasso.fit <-huge(cov(X), method = "glasso", nlambda=40,lambda.min.ratio = 0.01)
omega_path_glasso0 <- omega_path_glasso <- array(NA, dim = c(M, M, 40))
for(i in 1:dim(omega_path_glasso)[3]){
	omega_path_glasso[,,i] <- as.matrix(glasso.fit$icov[[i]])
	omega_path_glasso0[,,i] <- cov2cor(as.matrix(glasso.fit$icov[[i]]))
	diag(omega_path_glasso0[, , i]) <- NA
}
#########################################################################
pdf("../figures/illustration-emgs.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
omega_path0 <- omega_path
for(i in 1:dim(omega_path0)[1]){
	omega_path0[i, , ] <- cov2cor(omega_path0[i, , ])
	diag(omega_path0[i, , ]) <- NA
}
plot(v0s, rep(0, length(v0s)), ylim = range(omega_path0, na.rm = TRUE), col = "white", xlab = "v0", main = "EMGS", ylab = "Partial correlation")
abline(h = 0, col = "red", lty = 2)
abline(h = 0.5, col = "red", lty = 2)
for(i in 2:M){
	for(j in 1:(i-1)){
		if(Prec[i, j] != 0){
			lines(v0s, omega_path0[, i, j] * (prob_path[, i, j] > 0.5), type = "l", col = "darkblue")
		}
		if(Prec[i, j] == 0){
			lines(v0s, omega_path0[, i, j] * (prob_path[, i, j] > 0.5), type = "l", col = "darkgreen", lty = 4)
		}
	}
}
legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("Edges", "Non-edges"))

plot(glasso.fit$lambda, rep(0, length(glasso.fit$lambda)), ylim = range(omega_path_glasso0, na.rm = TRUE), col = "white", xlab = "rho", main = "Graphical Lasso", ylab = "Partial correlation")
abline(h = 0, col = "red", lty = 2)
abline(h = 0.5, col = "red", lty = 2)
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
legend("bottomright", lty = c(1, 4), col = c("darkblue", "darkgreen"), c("Edges", "Non-edges"))
dev.off()