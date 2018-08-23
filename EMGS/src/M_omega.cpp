#include "M_omega.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;

mat M_omega(int N, int M, mat &omega, mat &S, double lambda, mat &Ed, uvec pseq) {
	int i;
	double v;
	vec u(M-1);
	uvec remove_i(M-1);
	uvec left_i(1);
	mat out_mi_mi(M-1, M-1);
	mat Ed_i_mi(1, M-1);
	mat S_i_mi(M-1, 1);

	mat out = omega;
	for(i = 0; i < M; i++){
		remove_i = find(pseq != i);
		left_i(0) = i;
		out_mi_mi = out.submat(remove_i, remove_i);
		Ed_i_mi = Ed.submat(left_i, remove_i);
		S_i_mi = S.submat(remove_i, left_i);

		v = N / (lambda + S(i, i));
		mat invsub = inv(out_mi_mi);
		u = -inv((lambda + S(i, i)) * invsub + diagmat(Ed_i_mi)) * S_i_mi;

		out.submat(left_i, remove_i) = u.t();
		out.submat(remove_i, left_i) = u;
		mat uinvu = u.t() * invsub * u;
		out(i, i) = v + uinvu(0, 0);
	}
	// for(i = 0; i < M; i++){
	// 	out.swap_rows(i, M-1);
	// 	out.swap_cols(i, M-1);
	// 	Ed.swap_rows(i, M-1);
	// 	Ed.swap_cols(i, M-1);
	// 	S.swap_rows(i, M-1);
	// 	S.swap_cols(i, M-1);
		
	// 	v = N / (lambda + S(M-1, M-1));
	// 	mat invsub = inv(out.submat(0, 0, M-2, M-2));
	// 	u = -inv((lambda + S(M-1, M-1)) * invsub + diagmat(Ed.submat(M-1, 0, M-1, M-2))) * S.submat(0, M-1, M-2, M-1);

	// 	out.submat(M-1, 0, M-1, M-2) = u.t();
	// 	out.submat(0, M-1, M-2, M-1) = u;
	// 	mat uinvu = u.t() * invsub * u;
	// 	out(M-1, M-1) = v + uinvu(0, 0);

	// 	out.swap_rows(i, M-1);
	// 	out.swap_cols(i, M-1);
	// 	Ed.swap_rows(i, M-1);
	// 	Ed.swap_cols(i, M-1);
	// 	S.swap_rows(i, M-1);
	// 	S.swap_cols(i, M-1);
	// }

	return(out);
}

// for(i in 1:M){
// 	v <- N / (lambda + S[i, i])
// 	u <- -solve((lambda + S[i, i]) * solve(omega_new[-i, -i]) + diag(Ed[i, -i])) %*% S[i, -i]
// 	omega_new[i, -i] <- u 
// 	omega_new[-i, i] <- omega_new[i, -i]
// 	omega_new[i, i] <- v + t(u) %*% solve(omega_new[-i, -i]) %*% u
// }