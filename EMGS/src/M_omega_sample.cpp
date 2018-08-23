#include "M_omega_sample.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;

vec mvrnormArma(vec mu, mat sigma) {
   int ncols = sigma.n_cols;
   mat Y = arma::randn(1, ncols);
   mat out = Y * arma::chol(sigma);
   out.each_row() += mu.t();
   return out.row(0).t();
}

mat M_omega_sample(int N, int M, mat &omega, mat &S, double lambda, mat &Ed, uvec pseq) {

	int i;
	double v;
	vec u(M-1);
	vec u2(M-1);
	uvec remove_i(M-1);
	uvec left_i(1);
	mat out_mi_mi(M-1, M-1);
	mat Ed_i_mi(1, M-1);
	mat S_i_mi(M-1, 1);


	mat out = omega;
	mat Cmat(M-1, M-1);

	for(i = 0; i < M; i++){
		remove_i = find(pseq != i);
		left_i(0) = i;
		out_mi_mi = out.submat(remove_i, remove_i);
		Ed_i_mi = Ed.submat(left_i, remove_i);
		S_i_mi = S.submat(remove_i, left_i);
		
		v = N / (lambda + S(i, i));
		mat invsub = inv(out_mi_mi);
		Cmat = inv((lambda + S(i, i)) * invsub + diagmat(Ed_i_mi));
		u = -Cmat * S_i_mi;
		u = mvrnormArma(u, Cmat);
		
		out.submat(left_i, remove_i) = u.t();
		out.submat(remove_i, left_i) = u;
		mat uinvu = u.t() * invsub * u;
		out(i, i) = v + uinvu(0, 0);

	}
	return(out);
}

// for(i in 1:M){
// 	v <- N / (lambda + S[i, i])
// 	u <- -solve((lambda + S[i, i]) * solve(omega_new[-i, -i]) + diag(Ed[i, -i])) %*% S[i, -i]
// 	omega_new[i, -i] <- u 
// 	omega_new[-i, i] <- omega_new[i, -i]
// 	omega_new[i, i] <- v + t(u) %*% solve(omega_new[-i, -i]) %*% u
// }