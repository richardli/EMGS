
###########################################################################
##  Application: Breakfast at Frats
##  - Data available online at: 
############################################################################
## Preprocess data
library(EMGS)
library(plyr)
library(huge)
library(corrplot)
data <- read.csv("../data/+Dunnhumby/dunnhumby _Breakfast-at-the-Frat/sales.csv")
labels <- read.csv("../data/+Dunnhumby/dunnhumby _Breakfast-at-the-Frat/prod.csv")
weeks <- as.character(unique(data[, 1]))
weeks <- sort(as.Date(weeks, format='%d-%b-%y'))
prods <- unique(data[, 3])
sales0 <- ddply(data, c("WEEK_END_DATE", "UPC"), function(x) colSums(x[c("UNITS", "SPEND")]))
sales0[,1] <- as.Date(as.character(sales0[,1]), , format='%d-%b-%y')

sales <- matrix(1, length(weeks), length(prods))
colnames(sales) <- prods
rownames(sales) <- weeks
price <- sales
for(i in 1:dim(sales0)[1]){
	w <- which(weeks == sales0[i, 1])
	s <- which(prods == sales0[i, 2])
	sales[w, s] = sales0[i, 4]
	price[w, s] = sales0[i, 4] / sales0[i, 3]
}
cats <- labels[match(prods, labels[,1]), 4]

## Detrend sales
Y0<- as.vector(log(sales)[-1, ])
X <-  as.vector(log(sales)[-156, ])
res <- residuals(lm(Y0~X))
res <- matrix(res, 155, 55)

## model fitting
set.seed(1234)
Y0 <- Y <- scale(res)
Y[sample(1:length(Y0), 100)] <- NA

v0s <-seq(0.01, 1, len = 20)
v1 <- 100
a <- 1
b <- 10
lambda <- 1
groups <- as.numeric(cats)-1
outcv <- cv.EMGS((Y), v0s, v1, lambda, a, b, epsilon=0.001, K=4, maxitr = 1e4, copula = FALSE, group = NULL, copula_itr=0)
outcv1 <- cv.EMGS((Y), v0s, v1, lambda, a, b, epsilon=0.001, K=4, maxitr = 1e4, copula = FALSE, group = groups, copula_itr=0)
glasso.fit <-huge(Y, method = "glasso", nlambda=40, lambda.min.ratio = 0.01)
outcv2 <- cv.glasso(Y, rholist = glasso.fit$lambda, K=5)


## Visualization of correlation matrix
corrwraper <- function(mat, title, group, lim){
	if(!is.null(group)){
		order <- order(as.character(group))
	}else{
		order <- 1:dim(mat)[1]
	}
	corrplot2(mat[order, order],  tl.pos = 'n', title = title, mar = c(0,0,1,0), method="color", bg = "gray70", na.label.col = "gray70", is.corr=F, cl.lim = lim)
	if(!is.null(group)){
		corrRect(table(as.character(group)))
	}
}

## Visualization
lim <- range(c(outcv$fit.min$omega[,,1], 
			   outcv1$fit.min$omega[,,1], 
			   as.vector(outcv2$fit.min$icov[[1]])))
par(mfrow = c(3, 3), oma=c(0,0,2,0)+1)
prec1 <- (as.matrix(outcv$fit.min$omega[,,1]))
# diag(prec1) <- 1
corrwraper(prec1, title = "Exchangeable prior", group=cats, lim = lim)
prec1[outcv$fit.min$path[,,1] <= 0.5] <- 0
prec1[prec1 == 0] <- NA
# diag(prec1) <- 1
corrwraper(prec1, title = "Exchangeable prior", group=cats, lim = lim)

ll <- apply(outcv$ll, 1, mean)
plot(v0s, ll, xlab = "v0", ylab = "negative log likelihood", type = "b", main = "Exchangeable prior")
points(outcv$v0.min, ll[which(v0s == outcv$v0.min)], col = "red", pch = 20)
abline(v = outcv$v0.min, lty=3, col="grey20")


prec2<- (as.matrix(outcv1$fit.min$omega[,,1]))
# diag(prec2) <- 1
corrwraper(prec2, title = "Structured prior", group=cats, lim = lim)
prec2[outcv1$fit.min$path[,,1] <= 0.5] <- 0
prec2[prec2 == 0] <- NA
# diag(prec2) <- 1
corrwraper(prec2, title = "Structured prior", group=cats, lim = lim)

ll1 <- apply(outcv1$ll, 1, mean)
plot(v0s, ll1, xlab = "v0", ylab = "negative log likelihood", type = "b", main = "Structured prior")
points(outcv1$v0.min, ll1[which(v0s == outcv1$v0.min)], col = "red", pch = 20)
abline(v = outcv1$v0.min, lty=3, col="grey20")


prec3 <- (as.matrix(outcv2$fit.min$icov[[1]]))
# diag(prec3) <- 1
corrwraper(prec3, title = "Graphical lasso", group=cats, lim = lim)
prec3[prec3 == 0] <- NA
# diag(prec3) <- 1
corrwraper(prec3, title = "Graphical lasso", group=cats, lim = lim)

ll2 <- apply(outcv2$ll, 1, mean)
rholist <- outcv2$rholist
plot(log(rholist), ll2, xlab = "log lambda", ylab = "negative log likelihood", type = "b", main = "Graphical lasso")
points(log(outcv2$rho.min), ll2[which(rholist == outcv2$rho.min)], col = "red", pch = 20)
abline(v = log(outcv2$rho.min), lty=3, col="grey20")


# Missing data imputation Toy example
# set.seed(123)
# library(EMGS)
# library(huge)
# library(tmvtnorm)
# # sourceCpp("../EMGS/src/_EMGS.cpp")
# N <- 100
# M <- 10
# Prec <- matrix(0, M, M)
# diag(Prec) <- 1
# for(i in 2:M){
# 	Prec[i-1, i] <- 0.5
# 	Prec[i, i - 1] <- 0.5
# }
# nedge <- (sum(Prec != 0)-M)/2
# diag(Prec) <- 1
# # Prec <- Prec * 2
# Sigma <- solve(Prec)
# X0 <- rtmvnorm(n = N, mean = rep(0, M), sigma = Sigma)
# X <- X0
# X[sample(1:(N*M), 100)] <- NA
# ################################################################
# v0s <- seq(0.001, 0.5, len = 40)
# v1 <- 100
# epsilon <- 1e-6
# a <- 2
# b <- M
# lambda <- 1
# omega_path <- array(0, dim = c(length(v0s), M, M))
# prob_path <- array(0, dim = c(length(v0s), M, M))
# fit.cv <- cv.EMGS(X, v0s, v1, lambda, a, b, epsilon, copula=FALSE, copula_itr = 0)
# cols = matrix(1, N, M)
# cols[is.na(X)] <- 2
# plot(X0, fit.cv$fit.min$X, col = cols)



