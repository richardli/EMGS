#ifndef _EMGS_GET_MISSING_H
#define _EMGS_GET_MISSING_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat Get_missing(mat &X, mat &omega, mat &S);

#endif