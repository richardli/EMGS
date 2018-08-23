#include "_EMGS_rev.h"
#include "Get_expectation.h"
#include "M_omega_sample.h"
#include "M_omega.h"
#include "M_tau.h"
#include "Estimate_S.h"

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 


using namespace std;
using namespace Rcpp;
using namespace arma;

// debug: [[Rcpp::export]]
SEXP EMGS_rev(SEXP X_r, SEXP S_r, SEXP v0_r, SEXP v1_r, SEXP lambda_r, SEXP a_r, SEXP b_r, SEXP epsilon_r, SEXP verbose_r, SEXP maxitr_r, SEXP savepath_r, SEXP copula_r, SEXP Sitr_r, SEXP ranks_r, SEXP mranks_r, SEXP Z_init_r, SEXP thin_r, SEXP exist_group_r, SEXP group_r, SEXP a_tau_r, SEXP b_tau_r) {
    // organize input
    mat X = as<mat>(X_r);
    mat S = as<mat>(S_r);
    mat Scurrent = as<mat>(S_r);
    vec v0s = as<vec>(v0_r);
    double v1 = as<double>(v1_r);
    double lambda = as<double>(lambda_r);
    double a = as<double>(a_r);
    double b = as<double>(b_r);
    double epsilon = as<double>(epsilon_r);
    int maxitr = as<int>(maxitr_r);
    bool verbose = as<bool>(verbose_r);
    bool savepath = as<bool>(savepath_r);
    bool copula = as<bool>(copula_r);
    int Sitr = as<int>(Sitr_r);
    // int thin = as<int>(thin_r);
    int exist_group = as<int>(exist_group_r);
    uvec group = as<uvec>(group_r);;
    double a_tau = as<double>(a_tau_r);
    double b_tau = as<double>(b_tau_r);

    // take in the ranks of the variables, start from 0 being smallest
    umat ranks = as<umat>(ranks_r);
    uvec mranks = as<uvec>(mranks_r);
    mat ZZ = as<mat>(Z_init_r);

    // get dimensions
    const int L = v0s.n_elem;
	const int M = X.n_cols;
    const int N = X.n_rows;

    // get transformed parameters
    // mat S = X.t() * X;
	// mat Xt = X.t(); 

    // initialize output
    cube omega(M, M, L);
    cube EZ(M, M, L);
    cube Ed(M, M, L);
    cube tau_out(M, M, L);
    vec thetas(L);
 	vec itr(L);

    // temporary output to be removed later
    cube Edpath(M, M, 500);

    // initialize intermediate values
    int niter, i, count, tt;
    double v0, eps, eps0, EZsum, tempo;
    mat u(M-1, M-1);
	cube tmp(M, M, 2);
    uvec pseq(M);
    for(i = 0; i < M; i++){
     pseq(i) = i;
    }


    for(i = 0; i < L; i++){
    	v0 = v0s[i];
		mat v0eye = eye<mat>(M, M);
    	v0eye = v0eye * v0;
        mat omega_new = N * inv(S + v0eye);
        mat omega_update = omega_new;
        mat omega_current = omega_new;
        mat omega_2_current = omega_new;
        mat omega_2_update = omega_new;
        mat omega_mean_old = eye<mat>(M, M);
        mat omega_mean = eye<mat>(M, M);
	    mat Z = eye<mat>(M, M);
	   
        mat EZ_new = eye<mat>(M, M);
        mat EZ_update = eye<mat>(M, M);
        double theta = a / (a + b);
        eps = epsilon + 1;
        eps0 = epsilon + 1;
    	niter = 1;
      	mat Ed_new = eye<mat>(M, M);
        mat tau = ones<mat>(M,M);
		count = 0;

    	while((eps>epsilon | eps0 > epsilon) & (count < maxitr)){
            
            tempo = 1 / (count + 1.0);
 
            // SAE-step, calculating 
            //   1. sqrt(E(omega^2))
            //   2. E(omega)
            omega_update = omega_current * (1-tempo);
            omega_2_update = square(omega_2_current) * (1-tempo);
            for(tt=0; tt < 1; tt++){
                omega_new = M_omega_sample(N, M, omega_current, S, lambda, Ed_new, pseq);
                omega_update += omega_new * tempo / (1.0);
                omega_2_update += square(omega_new) * tempo / (1.0);
            }
            // tmp =  M_omega_sample(N, M, omega_current, S, lambda, Ed_new);
            // omega_update = tmp.slice(0);
            // omega_2_update = tmp.slice(1);

            omega_2_current = sqrt(omega_2_update);
            omega_current = omega_update;
            omega_new = omega_update;

            // M-step
			EZ_update = Get_expectation(M, omega_2_current, v0, v1, theta);
            EZ_update.elem( find(EZ_update > 0.5) ).ones();
            EZ_update.elem( find(EZ_update < 0.5) ).zeros();
			//Rcout << "EZ_update = \n" << EZ_update << endl;

			EZsum = accu(EZ_update) / 2; // diagonals are zero
			//Rcout << "EZsum = \n" << EZsum << endl;

			Ed_new = EZ_update / v1 / v1 + (1 - EZ_update) / v0 / v0;
            Ed_new = Ed_new % tau;

            if(copula){
               // this function also updates Z
               //Rcout << "omega_new = \n" << omega_new << endl;
               //Rcout << "ZZ0 = \n" << ZZ.row(0) << endl;
               S = Estimate_S(N, M, omega_2_current, Sitr, ranks, mranks, ZZ);
               S = Scurrent*(1-tempo) + S*tempo;
               Scurrent = S;
               //Rcout << "S new = \n" << S << endl;
            }

            // // still a M-step?
            // omega_update = M_omega(N, M, omega_new, S, lambda, Ed_new);
            // omega_current = omega_update;

			
            if(exist_group > 0){
                tau = M_tau(exist_group, M, group, Ed_new,omega_2_current, tau, a_tau, b_tau);
            }
        	// Rcout << "omega_update = \n" << omega_update << endl;

			theta = (a - 1 + EZsum) / (a + b + M * (M -1) / 2 - 2);

          

        	// Rcout << "theta = \n" << theta << endl;

            // if((count > thin) && copula){
            //     omega_mean = omega_new * (count - thin) + omega_mean;
            //     omega_mean /= (count - thin);
            //     eps = max(max(abs(omega_mean - omega_mean_old)));
            //     omega_mean_old = omega_mean;
            // }else{
            //     eps = max(max(abs(omega_update - omega_new)));
            // }
            // eps = max(max(abs(omega_update - omega_new)));
            eps0 = eps;
            eps = sum(sum(abs(EZ_new- EZ_update)));
            EZ_new = EZ_update;
        	count++;
        	if(verbose){
	        	Rcout << "Itr = " << count << " Max diff = " << eps << endl;
        	}else{
        		Rcout << ".";
        	}
            // // temporary output to be removed later
            // if(i == 0 & count < 500 & savepath){
            //     Edpath.slice(count) = Ed_new;
            // }
    	}
        // if(copula){
        //     omega.slice(i) = omega_mean;
        //     EZ.slice(i) = Get_expectation(M, omega_mean, v0, v1, theta);
        //     Ed.slice(i) = EZ.slice(i) / v1 / v1 + (1 - EZ_update) / v0 / v0;
        //     EZsum = accu(EZ.slice(i)) / 2; 
        //     theta = (a - 1 + EZsum) / (a + b + M * (M -1) / 2 - 2);
        // }else{
            omega.slice(i) = omega_new;
            EZ.slice(i) = EZ_update;
            Ed.slice(i) = Ed_new;
            tau_out.slice(i) = tau;
            thetas(i) = theta;
        // } 	
    	itr(i) = count;
    	Rcout << "v0 = " << v0 << " done" << endl;
    }
    
    List list;
    list["S"] = S;
    list["EZ"] = EZ;
    list["Ed"] = Ed;
    list["omega"] = omega;
	list["theta"] = thetas;    
	list["itr"] = itr;  
    list["tau"] = tau_out;
  
    if(savepath){
        list["Edpath"] = Edpath;
    }
    return list;
}

