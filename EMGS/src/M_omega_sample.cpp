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

mat M_omega_sample(int N, int M, mat &omega, mat &S, double lambda, mat &Ed) {
	int i;
	double v;
	vec u(M-1);
	vec u2(M-1);
	mat out = omega;
	mat Cmat(M-1, M-1);

	// mat out2 = omega;
	// cube out3(M, M, 2);

	for(i = 0; i < M; i++){
		out.swap_rows(i, M-1);
		out.swap_cols(i, M-1);
		Ed.swap_rows(i, M-1);
		Ed.swap_cols(i, M-1);
		S.swap_rows(i, M-1);
		S.swap_cols(i, M-1);
		
		v = N / (lambda + S(M-1, M-1));
		mat invsub = inv(out.submat(0, 0, M-2, M-2));
		Cmat = inv((lambda + S(M-1, M-1)) * invsub + diagmat(Ed.submat(M-1, 0, M-1, M-2)));
		u = -Cmat * S.submat(0, M-1, M-2, M-1);
		u = mvrnormArma(u, Cmat);
		
		// u2 = Cmat.diag() + square(u);
		// out2.swap_rows(i, M-1);
		// out2.swap_cols(i, M-1);
		// out2.submat(M-1, 0, M-1, M-2) = u2.t();
		// out2.submat(0, M-1, M-2, M-1) = u2;
		// out2.swap_rows(i, M-1);
		// out2.swap_cols(i, M-1);

		out.submat(M-1, 0, M-1, M-2) = u.t();
		out.submat(0, M-1, M-2, M-1) = u;
		mat uinvu = u.t() * invsub * u;
		out(M-1, M-1) = v + uinvu(0, 0);

		out.swap_rows(i, M-1);
		out.swap_cols(i, M-1);
		Ed.swap_rows(i, M-1);
		Ed.swap_cols(i, M-1);
		S.swap_rows(i, M-1);
		S.swap_cols(i, M-1);
	}
	// out3.slice(0) = out;
	// out3.slice(1) = out2;

	return(out);
}

// for(i in 1:M){
// 	v <- N / (lambda + S[i, i])
// 	u <- -solve((lambda + S[i, i]) * solve(omega_new[-i, -i]) + diag(Ed[i, -i])) %*% S[i, -i]
// 	omega_new[i, -i] <- u 
// 	omega_new[-i, i] <- omega_new[i, -i]
// 	omega_new[i, i] <- v + t(u) %*% solve(omega_new[-i, -i]) %*% u
// }