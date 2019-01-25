#include "Get_missing.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;


// debug[[Rcpp::export]]

mat Get_missing(mat &X, mat &omega, mat &S){
	int i;
	uvec ii(1);

	for(i = 0; i < X.n_rows; i++){
		uvec miss = arma::find_nonfinite(X.row(i));
		if(miss.n_elem > 0){
			ii(0) = i;
			uvec remain = arma::find_finite(X.row(i));
			mat tmp = inv(omega.submat(miss, miss));
			X.submat(ii, miss) = -(tmp * omega.submat(miss, remain) * X.submat(ii, remain).t()).t();
			S.submat(miss, miss) = S.submat(miss, miss)	+ tmp;
		}
	}
	return(X);
}



// sourceCpp("src/Get_missing.cpp")