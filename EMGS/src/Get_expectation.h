#ifndef _EMGS_Get_expectation_H
#define _EMGS_Get_expectation_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat Get_expectation(int M, mat &omega, double v0, double v1, double theta) ;

#endif