#include "Get_expectation.h"
#include "density_norm.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;


mat Get_expectation(int M, mat &omega, double v0, double v1, double theta) {
	mat aa = zeros<mat>(M, 2);
	mat EZ = zeros<mat>(M, M);
	int i;
	vec abmax(M);
	
	for(i = 0; i < M; i++){
			vec tmp = omega.col(i);
			aa.col(1) = density_norm_log(tmp, 0.0, v1) + log(theta);
			aa.col(0) = density_norm_log(tmp, 0.0, v0) + log(1 - theta);
			abmax = max(aa, 1); // max for each row
			// Rcout << "aa inside = \n" << aa << endl;

			EZ.col(i) = exp(aa.col(1) - abmax) / (exp(aa.col(0) - abmax) + exp(aa.col(1) - abmax));
			EZ.row(i) = EZ.col(i).t();
			EZ(i, i) = 0;
	}

    return EZ;
}

 