#ifndef _EMGS_M_tau_H
#define _EMGS_M_tau_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat M_tau(int exist_group, int M, uvec &group, mat &Ed_new, mat &omega_update, mat &tau, double a_tau, double b_tau);

#endif