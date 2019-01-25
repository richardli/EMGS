#ifndef _EMGS_DNORM_H
#define _EMGS_DNORM_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
vec density_norm(vec &x, double mu, double sigma);
vec density_norm_log(vec &x, double mu, double sigma);
vec density_norm_log_var(vec &x, double mu, vec &sigma2);
vec pnorm_log(vec x);
vec pnorm(vec x);
double erfc_log(double x);
#endif
