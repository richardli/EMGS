#ifndef _EMGS_M_omega_sample_H
#define _EMGS_M_omega_sample_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat M_omega_sample(int N, int M, mat &omega, mat &S, double lambda, mat &Ed) ;

#endif