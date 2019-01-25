###################################################################
## Gaussian case
load(paste("../data/simG-results/p", p, "/graphfit", n, p, graph, itr, "gaussian.rda", sep="_"))
Prec <- fit$Prec
G <- fit$G
fit.emgs <- fit$cv.out$fit.min
fit.glasso <- fit$cv.out.npn$fit.min
fit.glasso.ric <- fit$out.ric
fit.glasso.star <- fit$out.star
fit.glasso.xz <- NULL

prec.emgs <- fit.emgs$omega[,,1]
prec.glasso <- as.matrix(fit.glasso$icov[[1]])
prec.glasso.ric <- as.matrix(fit.glasso.ric$opt.icov)
prec.glasso.star <- as.matrix(fit.glasso.star$opt.icov)
prec.glasso.xz <- NULL

norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star), Prec)
roc.emgs <- getROCthres(prec.emgs, G)
roc.glasso <- getROCthres(prec.glasso, G)
roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
roc.glasso.star <- getROCthres(prec.glasso.star, G)
auc.emgs <- roc.emgs$AUC
auc.glasso <- roc.glasso$AUC
auc.glasso.ric <- roc.glasso.ric$AUC
auc.glasso.star <- roc.glasso.star$AUC
graph.emgs <- Threshold(prec.emgs, sum(G)/2)
graph.glasso <- Threshold(prec.glasso, sum(G)/2)
graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
f1 <- getF1(list(graph.emgs, graph.glasso, graph.glasso.ric, graph.glasso.star), G)
  # AUC
AUC_lam <- c(getROCpath(fit$out$path, G)$AUC, getROCpath(fit$out.npn$path, G)$AUC, NA, NA)
# F1
graph.raw <- fit$cv.out$fit.min$path[,,1]
graph.gl  <- prec.glasso
graph.gl[graph.gl != 0] <- 1        
graph.ric  <- prec.glasso.ric
graph.ric[graph.ric != 0] <- 1 
graph.star <- prec.glasso.star
graph.star[graph.star != 0] <- 1 
F1_raw <-  getF1(list(graph.raw, graph.gl, graph.ric, graph.star), G)
eval.g <- list(norms=norms, auc = c(auc.emgs, auc.glasso, auc.glasso.ric, auc.glasso.star), f1 = f1, Auc_lam = AUC_lam, f1_raw = F1_raw)



###################################################################
## Non-Gaussian case
load(paste("../data/simM-results/p", p, "/graphfit", n, p, graph, itr, "mixed.rda", sep="_"))
Prec <- fit$Prec
G <- fit$G
fit.emgs <- fit$cv.out$fit.min
fit.glasso <- fit$cv.out.npn$fit.min
fit.glasso.ric <- fit$out.ric
fit.glasso.star <- fit$out.star
fit.glasso.xz <- fit$cv.out.xz$fit.min

prec.emgs <- fit.emgs$omega[,,1]
prec.glasso <- as.matrix(fit.glasso$icov[[1]])
prec.glasso.ric <- as.matrix(fit.glasso.ric$opt.icov)
prec.glasso.star <- as.matrix(fit.glasso.star$opt.icov)
prec.glasso.xz <- as.matrix(fit.glasso.xz$icov[[1]])

norms <- getNorms(list(prec.emgs, prec.glasso, prec.glasso.ric, prec.glasso.star, prec.glasso.xz), Prec)
roc.emgs <- getROCthres(prec.emgs, G)
roc.glasso <- getROCthres(prec.glasso, G)
roc.glasso.ric <- getROCthres(prec.glasso.ric, G)
roc.glasso.star <- getROCthres(prec.glasso.star, G)
roc.glasso.xz <- getROCthres(prec.glasso.xz, G)
auc.emgs <- roc.emgs$AUC
auc.glasso <- roc.glasso$AUC
auc.glasso.ric <- roc.glasso.ric$AUC
auc.glasso.star <- roc.glasso.star$AUC
auc.glasso.xz <- roc.glasso.xz$AUC
graph.emgs <- Threshold(prec.emgs, sum(G)/2)
graph.glasso <- Threshold(prec.glasso, sum(G)/2)
graph.glasso.ric <- Threshold(prec.glasso.ric, sum(G)/2)
graph.glasso.star <- Threshold(prec.glasso.star, sum(G)/2)
graph.glasso.xz <- Threshold(prec.glasso.xz, sum(G)/2)
f1 <- getF1(list(graph.emgs, graph.glasso, graph.glasso.ric, graph.glasso.star, graph.glasso.xz), G)
# AUC
AUC_lam <- c(getROCpath(fit$out$path, G)$AUC, getROCpath(fit$out.npn$path, G)$AUC, NA, NA, getROCpath(fit$out.xz$path, G)$AUC)
# F1
graph.raw <- fit$cv.out$fit.min$path[,,1]
graph.gl  <- prec.glasso
graph.gl[graph.gl != 0] <- 1        
graph.ric  <- prec.glasso.ric
graph.ric[graph.ric != 0] <- 1 
graph.star <- prec.glasso.star
graph.star[graph.star != 0] <- 1 
graph.xz <- prec.glasso.xz
graph.xz[graph.xz != 0] <- 1 
F1_raw <-  getF1(list(graph.raw, graph.gl, graph.ric, graph.star, graph.xz), G)
eval.m <- list(norms=norms, auc = c(auc.emgs, auc.glasso, auc.glasso.ric, auc.glasso.star, auc.glasso.xz), f1 = f1, Auc_lam = AUC_lam, f1_raw = F1_raw)

out <- list(Gaussian = eval.g, Mixed = eval.m)
save(out, file = paste("../data/metrics/", n, p, graph, itr, "metrics.rda", sep="_"))     
print(out)