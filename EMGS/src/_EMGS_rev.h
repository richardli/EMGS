#ifndef _EMGS_EMGS_reverse_H
#define _EMGS_EMGS_reverse_H

#include <RcppArmadillo.h>

RcppExport SEXP EMGS_rev(SEXP X_r, SEXP S_r, SEXP v0_r, SEXP v1_r, SEXP lambda_r, SEXP a_r, SEXP b_r, SEXP epsilon_r, SEXP verbose_r, SEXP maxitr_r, SEXP savepath_r, SEXP copula_r, SEXP Sitr_r, SEXP ranks_r, SEXP mranks_r, SEXP Z_init_r, SEXP thin_r, SEXP exist_group_r, SEXP group_r, SEXP a_tau_r, SEXP b_tau_r);

#endif