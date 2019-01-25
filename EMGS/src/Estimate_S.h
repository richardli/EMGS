#ifndef _EMGS_Estimate_S_H
#define _EMGS_Estimate_S_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat Estimate_S(int N, int M, mat &omega, int Sitr, umat &ranks, uvec &mranks, mat &ZZ, uvec pseq) ;

#endif