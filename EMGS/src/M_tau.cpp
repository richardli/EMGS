#include "M_tau.h"
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;

mat M_tau(int exist_group, int M, uvec &group, mat &Ed_new, mat &omega_update, mat &tau, double a_tau, double b_tau) {
	int i, j;
	vec u(M-1);
	mat tau_compact= ones<mat>(exist_group, exist_group);
	mat size = eye<mat>(exist_group, exist_group);
	tau_compact = tau_compact * b_tau;

	for(i=0; i < M; i++){
		for(j=0; j < M; j++){
			if(i < j){
				tau_compact(group(i), group(j)) += 0.5 * omega_update(i, j) * omega_update(i, j) * Ed_new(i, j);
				size(group(i), group(j)) += 0.5;
				if(group(i) != group(j)){
					tau_compact(group(j), group(i)) += 0.5 * omega_update(j, i) * omega_update(j, i) * Ed_new(j, i);
					size(group(j), group(i)) += 0.5;					
				}
			}
		}
	}
	for(i=0; i < exist_group; i++){
		for(j=0; j < exist_group; j++){
			tau_compact(i,j) = (a_tau - 1 + size(i,j)) / tau_compact(i,j);
		}
	}
    // Rcout << "tau = " << tau_compact << endl;
    // Rcout << "size = " << size << endl;


	for(i=0; i < M; i++){
		for(j=0; j < M; j++){
			tau(i, j) = tau_compact(group(i), group(j)); 
		}
	}

	return(tau);
}
