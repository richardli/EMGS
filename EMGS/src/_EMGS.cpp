#include "_EMGS.h"
#include "Get_expectation.h"
#include "Get_missing.h"
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
SEXP _EMGS(SEXP X_r, SEXP S_r, SEXP v0_r, SEXP v1_r, SEXP lambda_r, SEXP a_r, SEXP b_r, SEXP epsilon_r, SEXP verbose_r, SEXP maxitr_r, SEXP savepath_r, SEXP copula_r, SEXP Sitr_r, SEXP ranks_r, SEXP mranks_r, SEXP Z_init_r, SEXP thin_r, SEXP exist_group_r, SEXP group_r, SEXP a_tau_r, SEXP b_tau_r, SEXP warm_r, SEXP Xt_r, SEXP St_r, SEXP ranks_t_r, SEXP mranks_t_r, SEXP Zt_init_r, SEXP missing_r) {
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
    bool warm = as<bool>(warm_r);
    bool has_missing = as<bool>(missing_r);
    mat X0 = as<mat>(X_r); // to hold the X with missing values in it
    mat S0 = as<mat>(S_r); // to hold the S with missing values = 0

    // take in the ranks of the variables, start from 0 being smallest
    umat ranks = as<umat>(ranks_r);
    uvec mranks = as<uvec>(mranks_r);
    mat ZZ = as<mat>(Z_init_r);

    // get dimensions
    const int L = v0s.n_elem;
	const int M = X.n_cols;
    const int N = X.n_rows;

    // testing data
    mat Xt = as<mat>(Xt_r);
    mat Xt0 = as<mat>(Xt_r);
    bool Xt_exist = (Xt.n_cols == X.n_cols);
    mat St = as<mat>(St_r);
    mat St0 = as<mat>(St_r);
    mat Stcurrent = as<mat>(St_r);
    umat ranks_t = as<umat>(ranks_t_r);
    uvec mranks_t = as<uvec>(mranks_t_r);
    mat ZZt = as<mat>(Zt_init_r);
    const int Nt = Xt.n_rows;
    cube Stout(M, M, L);

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
    int niter, i, count, ii, jj;
    double v0, eps, EZsum, tempo;
    // double v;
    mat u(M-1, M-1);
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
        if(warm & (i > 0)){
            for(ii=0; ii < M; ii++){
                for(jj=0; jj < M; jj++){
                    if(ii != jj & EZ(ii, jj, i-1) < 0.5){
                        omega_new(ii, jj) = 0;
                    }else{
                        omega_new(ii, jj) = omega(ii, jj, i-1);
                    }                     
                }
            }
            mat omega_update = omega_new;        
        }
        mat omega_mean_old = eye<mat>(M, M);
        mat omega_mean = eye<mat>(M, M);
	    mat Z = eye<mat>(M, M);
	    double theta = 0.5;
		eps = epsilon + 1;
    	niter = 1;
		mat EZ_new = eye<mat>(M, M);
		mat Ed_new = eye<mat>(M, M);
        mat tau = ones<mat>(M,M);
		count = 0;

        // re-initiate copula latent variables
        ZZ = as<mat>(Z_init_r);
        Scurrent = as<mat>(S_r);
        ZZt = as<mat>(Zt_init_r);
        Stcurrent = as<mat>(St_r);            
        

    	while(eps>epsilon & count < maxitr){
            
            tempo = 1 / (count + 1.0);

            if(has_missing){
               S = S0 * 1; 
               X = X0 * 1;
               // this function also updates S
               X = Get_missing(X, omega_new, S);
               if(Xt_exist){
                    Xt = Xt0 * 1;
                    St = St0 * 1;
                    Xt = Get_missing(Xt, omega_new, St0);
               }
            }
			// E-step
			EZ_new = Get_expectation(M, omega_new, v0, v1, theta, tau);
			EZsum = accu(EZ_new) / 2; // diagonals are zero

			Ed_new = EZ_new / v1 / v1 + (1 - EZ_new) / v0 / v0;
            Ed_new = Ed_new % tau;

            if(copula){
               // this function also updates Z
               S = Estimate_S(N, M, omega_new, Sitr, ranks, mranks, ZZ, pseq);
               S = Scurrent*(1-tempo) + S*tempo;
               Scurrent = S;
               if(Xt_exist){
                   St = Estimate_S(Nt, M, omega_new, Sitr, ranks_t, mranks_t, ZZt, pseq);
                   St = Stcurrent*(1-tempo) + St*tempo;
                   Stcurrent = St;                
               }
            }

			// M-step
			omega_update = M_omega(N, M, omega_new, S, lambda, Ed_new, pseq);
            if(exist_group > 0){
                tau = M_tau(exist_group, M, group, Ed_new,omega_update, tau, a_tau, b_tau);
            }

			theta = (a - 1 + EZsum) / (a + b + M * (M -1) / 2 - 2);

            // if((count > thin) && copula){
            //     omega_mean = omega_new * (count - thin) + omega_mean;
            //     omega_mean /= (count - thin);
            //     eps = max(max(abs(omega_mean - omega_mean_old)));
            //     omega_mean_old = omega_mean;
            // }else{
            //     eps = max(max(abs(omega_update - omega_new)));
            // }
            eps = max(max(abs(omega_update - omega_new)));
        	omega_new = omega_update;
        	count++;
        	if(verbose){
	        	Rcout << "Itr = " << count << " Max diff = " << eps << endl;
        	}else{
        		Rcout << ".";
        	}
    	}
            omega.slice(i) = omega_new;
            EZ.slice(i) = EZ_new;
            Ed.slice(i) = Ed_new;
            tau_out.slice(i) = tau;
            thetas(i) = theta;
            if(Xt_exist) Stout.slice(i) = St;
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
    list["St"] = Stout;
    list["X"] = X;
  
    if(savepath){
        list["Edpath"] = Edpath;
    }
    return list;
}

