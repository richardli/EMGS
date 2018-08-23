#ifndef _EMGS_M_omega_H
#define _EMGS_M_omega_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat M_omega(int N, int M, mat &omega, mat &S, double lambda, mat &Ed, uvec pseq) ;

#endif