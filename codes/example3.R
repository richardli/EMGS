
library(EMGS)
library(plyr)
library(huge)
library(corrplot)
library(reshape2)
source("corrplot2.R")

dat = read.csv("../data/Burke_Gilman_Trail_north_of_NE_70th_St_Bike_and_Ped_Counter.csv", check.names=F)
dat <- data.frame(dat, check.names = F)
dat$Date <- as.character(dat$Date)
dat$Date <-strptime(dat$Date, format = "%m/%d/%y %H:%M")
dat <- dat[dat$Date < strptime("1/1/15 0:00", format = "%m/%d/%y %H:%M"), ]
dat <- dat[!is.na(dat$Date), ]
dat$Day <- as.numeric(trunc(difftime(dat$Date, min(dat$Date), units = "day"))) + 1
dat$hour <- as.character(format(dat$Date,   "%H"))
p <- c("small", "large")[2]

if(p =="large"){
	##########################################################
	dat2 <- matrix(0, 365, 24*4)
	for(j in 3:6){
		tmp <- dcast(dat, Day~hour, value.var = colnames(dat)[j], drop = FALSE, fun.aggregate = function(x){x[1]})
		colnames(tmp)[-1] <- paste(colnames(dat)[j], 0:23, sep = "_")
		dat2[, (j-3)*24 + 1:24] <- as.matrix(tmp[, -1])
	}
}else{
	##########################################################
	dat2 <- matrix(0, 365, 24*2)
	dat$Ped = dat[,3] + dat[,4]
	dat$Bike = dat[,5] + dat[,6]
	for(j in 9:10){
		tmp <- dcast(dat, Day~hour, value.var = colnames(dat)[j], drop = FALSE, fun.aggregate = function(x){x[1]})
		colnames(tmp)[-1] <- paste(colnames(dat)[j], 0:23, sep = "_")
		dat2[, (j-9)*24 + 1:24] <- as.matrix(tmp[, -1])
	}
	##########################################################
}
dat2[is.na(dat2)] <- 0
groups <- rep(0:(dim(dat2)[2]/24-1), each = 24)
times <- rep(0:23, (dim(dat2)[2]/24))

####################################################
## Fit Model to the entire dataset
####################################################
Y <- log(dat2+1)
mean <- apply(Y, 2, mean, na.rm = TRUE)
for(i in 1:dim(Y)[2]){
	Y[, i] <- Y[, i] - mean[i]
}

v0s <-seq(0.01, 1, len = 20)
v1 <- 100
a <- 1
b <- 10
lambda <- 1
outcv <- cv.EMGS((Y), v0s, v1, lambda, a, b, epsilon=0.001, K=5, maxitr = 1e4, copula = FALSE, group = NULL, copula_itr=0)
outcv1 <- cv.EMGS((Y), v0s, v1, lambda, a, b, epsilon=0.001, K=5, maxitr = 1e4, copula = FALSE, group = groups, copula_itr=0)
glasso.fit <-huge(Y, method = "glasso", nlambda=20, lambda.min.ratio = 0.001)
outcv2 <- cv.glasso(Y, rholist = glasso.fit$lambda, K=5)


####################################################
## Visualization
####################################################

pdf("../figures/burke.pdf", width = 9, height = 6)
par(mfrow = c(2, 3), oma=c(0,0,1,0))
cats <- groups
prec1 <- (as.matrix(outcv$fit.min$omega[,,1]))
prec2<- (as.matrix(outcv1$fit.min$omega[,,1]))
prec3 <- (as.matrix(outcv2$fit.min$icov[[1]]))


diag(prec1) <- diag(prec2) <- diag(prec3) <- NA
lim <- range(c(prec1, prec2, prec3), na.rm = TRUE)
prec1[outcv$fit.min$path[,,1] <= 0.5] <- 0
prec1[prec1 == 0] <- NA
corrwraper(prec1, title =expression(paste("Exchangeable prior: ", hat(Omega))), group=cats, lim = lim)
prec2[outcv1$fit.min$path[,,1] <= 0.5] <- 0
prec2[prec2 == 0] <- NA
corrwraper(prec2, title = expression(paste("Structured prior: ", hat(Omega))), group=cats, lim = lim)
prec3[prec3 == 0] <- NA
corrwraper(prec3, title = expression(paste("Graphical lasso: ", hat(Omega))), group=cats, lim = lim)


lim <- range(c(solve(outcv$fit.min$omega[,,1]), 
			   solve(outcv1$fit.min$omega[,,1]), 
			   as.vector(solve(outcv2$fit.min$icov[[1]]))))
prec1 <- (as.matrix(outcv$fit.min$omega[,,1]))
corrwraper(solve(prec1), title = expression(paste("Exchangeable prior: ", hat(Omega)^-1)), group=cats, lim = lim)
prec2<- (as.matrix(outcv1$fit.min$omega[,,1]))
corrwraper(solve(prec2), title =expression(paste("Structured prior: ", hat(Omega)^-1)), group=cats, lim = lim)
prec3 <- (as.matrix(outcv2$fit.min$icov[[1]]))
corrwraper(solve(prec3), title = expression(paste("Graphical lasso: ", hat(Omega)^-1)), group=cats, lim = lim)


dev.off()


####################################################
## Missing Data
####################################################
set.seed(1234)
Nitr <- 100
mse <- matrix(NA, Nitr, 4)
for(itr in 1:Nitr){
	Ymiss <- log(dat2+1)
	sample <- TRUE
	while(sample){
		miss <- sample(1:dim(dat2)[1], dim(dat2)[1]/2)
		mi <- sample(1:dim(Ymiss)[2], dim(Ymiss)[2]/2)
		sample <-  det(cov(Ymiss[-miss, ])) <= 0
	}
	ob <- (1:(dim(Ymiss)[2]/2))[-mi]
	X0 <- Ymiss[miss, mi]
	Ymiss[miss, mi] <- NA
	mean <- apply(Ymiss, 2, mean, na.rm = TRUE)
	for(i in 1:dim(Ymiss)[2]){
		Ymiss[, i] <- Ymiss[, i] - mean[i]
		if(i %in% mi){
			X0[, which(mi == i)] <- X0[, which(mi == i)] - mean[i]
		}
	}
	glasso.fit <-huge(cov(Ymiss[-miss, ]), method = "glasso", nlambda=20, lambda.min.ratio = 0.001)
	outcv2 <- cv.glasso(Ymiss[-miss, ], rholist = glasso.fit$lambda, K=5)

	outcv <- cv.EMGS((Ymiss), v0s, v1, lambda, a, b, epsilon=0.001, K=5, maxitr = 1e4, copula = FALSE, group = NULL, copula_itr=0)
	outcv1 <- cv.EMGS((Ymiss), v0s, v1, lambda, a, b, epsilon=0.001, K=5, maxitr = 1e4, copula = FALSE, group = groups, copula_itr=0)

	c(which.min(apply(outcv$ll, 1, mean)), which.min(apply(outcv1$ll, 1, mean)), which.min(apply(outcv2$ll, 1, mean)))

	X1 <- outcv$fit.min$X[miss, mi]
	X2 <- outcv1$fit.min$X[miss, mi]
	tmp <- outcv2$fit.min$icov[[1]]
	X3 <- t(-solve(tmp[mi, mi]) %*% tmp[mi, ob] %*% t(Ymiss[miss, ob]))
	tmp <- solve(cov(Ymiss[-miss, ]))
	X4 <- t(-solve(tmp[mi, mi]) %*% tmp[mi, ob] %*% t(Ymiss[miss, ob]))
	mse[itr, ] <- c(mean((X0 - X1)^2), mean((X0 - X2)^2), mean((X0 - X3)^2),  mean((X0 - X4)^2))
	message(apply(mse, 2, mean, na.rm = TRUE))
	message(apply(mse, 2, sd, na.rm = TRUE))
}







# ll <- apply(outcv$ll, 1, mean)
# plot(v0s, ll, xlab = "v0", ylab = "negative log likelihood", type = "b", main = "Exchangeable prior")
# points(outcv$v0.min, ll[which(v0s == outcv$v0.min)], col = "red", pch = 20)
# abline(v = outcv$v0.min, lty=3, col="grey20")

# ll1 <- apply(outcv1$ll, 1, mean)
# plot(v0s, ll1, xlab = "v0", ylab = "negative log likelihood", type = "b", main = "Structured prior")
# points(outcv1$v0.min, ll1[which(v0s == outcv1$v0.min)], col = "red", pch = 20)
# abline(v = outcv1$v0.min, lty=3, col="grey20")

# ll2 <- apply(outcv2$ll, 1, mean)
# rholist <- outcv2$rholist
# plot(log(rholist), ll2, xlab = "log lambda", ylab = "negative log likelihood", type = "b", main = "Graphical lasso")
# points(log(outcv2$rho.min), ll2[which(rholist == outcv2$rho.min)], col = "red", pch = 20)
# abline(v = log(outcv2$rho.min), lty=3, col="grey20")
# which.min(ll2)

# t1 = as.matrix(outcv$fit.min$omega[,,1])
# t2 = as.matrix(outcv1$fit.min$omega[,,1])
# e1 = as.matrix(outcv$fit.min$path[,,1])
# e2 = as.matrix(outcv1$fit.min$path[,,1])
# diag(t1) <- diag(t2) <- diag(e1) <- diag(e2) <- NA
# outcv1$fit.min$tau_compact
# tt <- NULL
# for(i in 1:dim(e1)[1]){
# 	for(j in 1:dim(e1)[2]){
# 		if(i<j && e1[i, j] != e2[i, j]){
# 			tt <- rbind(tt, c(i, j, ifelse(e1[i, j] > e2[i, j], 1, 2)))
# 		}
# 	}
# }
# tt <- data.frame(tt)
# tt$g1 <- groups[tt[, 1]]
# tt$g2 <- groups[tt[, 2]]
# tt$g <- paste(tt$g1, tt$g2, sep = "-")
# table(tt$g[tt$X3 == 1])
# table(tt$g[tt$X3 == 2])

