#include "Estimate_S.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;

// debug [[Rcpp::export]]
mat Estimate_S(int N, int M, mat &omega, int Sitr, umat &ranks, uvec &mranks, mat &ZZ, uvec pseq){
	// cout << N << endl;
	// cout << M << endl;
	// cout << X.row(0) << endl;
	// cout << omega << endl;
	// cout << Sitr << endl;
	// cout << ranks.row(0) << endl;
	// cout << mranks << endl;
 //    cout << ZZ.row(0) << endl;

	int i, m, j, r, k, jj;
	mat S = zeros<mat>(M, M);
	//vec mu = zeros<vec>(M);
	cube Sinv = zeros<cube>(M-1, M-1, M);
	mat sigma = inv(omega); 
	vec diag = zeros<vec>(M);
	for(i = 0; i < M; i++){
		diag(i) = sigma(i,i);
	}
	for(i = 0; i < M; i++){
		for(j = 0; j < M; j++){
			//sigma(i, j) = sigma(i,j) / sqrt(diag(i) * diag(j));
		}
	}

	// intermediate values
	double tmpF;
	double lower, upper;
	mat mean = zeros<mat>(1,1);
	mat var = zeros<mat>(1,1);
	rowvec omega_j = zeros<rowvec>(M);
	vec omega_j_mj = zeros<vec>(M-1);
	vec inv_times_sigma_j_mj = zeros<vec>(M-1);
	rowvec ZZ_i = zeros<rowvec>(M);
	vec ZZ_i_mj = zeros<vec>(M-1);
	uvec mj = zeros<uvec>(M);
	mat ZZmean = zeros<mat>(N, M);


	// // calculate Sinv first
	// for(j = 0; j < M; j++){
	// 	mj = find(pseq != j);
	// 	Sinv.slice(j) = inv(sigma.submat(mj, mj));
	// 	//Rcout << Sinv.slice(j) << endl;
	// 	// omega.submat(0, 0, M-2, M-2) - omega.submat(0, M-1, M-2, M-1) * omega.submat(M-1, 0, M-1, M-2) / omega(M-1, M-1);
	// }
	

	for(m = 0; m < Sitr; m++){
	
        // vec indices = shuffle(pseq);

		for(jj = 0; jj < M; jj++){
			j = jj;

			vec tmp = zeros<vec>(mranks(j));
			for(r = 0; r < mranks(j); r++){
				lower = R::qnorm(1 / (N + 1.0), 0.0, 1.0, true, false);
				upper = R::qnorm(1 - 1 / (N + 1.0), 0.0, 1.0, true, false);
				// Find positions of the below and upper elements
				arma::uvec pos_lower = find(ranks.col(j) == r-1);
				arma::uvec pos_current = find(ranks.col(j) == r);
				arma::uvec pos_upper = find(ranks.col(j) == r+1);
				// Find range
				for(k = 0; k < pos_lower.n_elem; k++){
					lower = max(lower, ZZ(pos_lower(k), j));
				}
				for(k = 0; k < pos_upper.n_elem; k++){
					upper = min(upper, ZZ(pos_upper(k), j));
				}
				// Calculate values for conditional mean and SD
				omega_j = omega.row(j);
				mj = find(pseq != j);
				omega_j_mj = omega_j.elem(mj);
				// var = sigma(j, j) - sigma_j_mj.t() * Sinv.slice(j) * sigma_j_mj;
				var = 1/omega(j, j);
				inv_times_sigma_j_mj = -1/omega(j,j) * omega_j_mj;

				// Sample all variables
				for(k = 0; k < pos_current.n_elem; k++){
					ZZ_i = ZZ.row(pos_current(k));
					ZZ_i_mj = ZZ_i.elem(mj);	
					mean = ZZ_i_mj.t() * inv_times_sigma_j_mj;
					
					// Rcout << pos_current(k) << j << mean << "-" << var << endl;
					tmpF = runif(1, 
								R::pnorm(lower, mean(0,0), sqrt(var(0,0)), true, false),
								R::pnorm(upper, mean(0,0), sqrt(var(0,0)), true, false))(0);
					if(tmpF > 1-1e-10){
						ZZ(pos_current(k), j) = lower + (upper - lower)/4;
					}else if(tmpF < 1e-10){
						ZZ(pos_current(k), j) = upper - (upper - lower)/4;
					}else{				
						ZZ(pos_current(k), j) = R::qnorm(tmpF, mean(0,0), sqrt(var(0,0)), true, false);			
					}
					double tmp  =  ZZ(pos_current(k), j);
					if(R_IsNaN(tmp)){
						Rcout << ZZ(pos_current(k), j) <<"_"<< tmpF <<"_"<<  mean << "_" << var << "_" << lower << "_" << upper << endl;
						return(var);
					}		
				} 

			}
		}
		// Calculate Y^T * Y
		S += ZZ.t() * ZZ ;
		// Calculate ZZ mean
		// ZZmean += ZZ;
	}
	// ZZmean /= (Sitr + 0.0);
	// S = ZZmean.t() * ZZmean;
	S /= Sitr;
	//Rcout << S << endl;
	return(S);
}