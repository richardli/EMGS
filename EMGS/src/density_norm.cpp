//density_norm.cpp
#include "density_norm.h"
#include <math.h> 

vec density_norm(vec &x, double mu, double sigma){
  vec dens = 1/(sigma*sqrt(2*PI)) * ones<vec>(x.n_elem);    
  dens %= exp(-square(x-mu)/(2 * pow(sigma, 2))); 
  return dens;
}


vec density_norm_log(vec &x, double mu, double sigma){
  vec dens = -0.5*log(sigma*sigma*2*PI) * ones<vec>(x.n_elem);    
  dens += -square(x-mu)/(2 * pow(sigma, 2)); 
  return dens;
}

vec pnorm(vec x){

  vec value=x/(sqrt(2));
  for (int i=0;i< x.n_elem;i++){
    value[i]=erf(value[i]);
    value[i]*=0.5;
    value[i]+=0.5;
  }
  return value;
}



vec pnorm_log(vec x){

  vec value=-x/(sqrt(2));
  for (int i=0;i< x.n_elem;i++){
    // value[i]=0.5+0.5*erf(value[i]);
    //value[i]=0.5*erfc(value[i]);
    //value[i]=log(value[i]);
    value[i]=log(0.5)+erfc_log(value[i]);
    
}
  return value;
}


double erfc_log(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y =  log (((((a5*t + a4)*t) + a3)*t + a2)*t + a1) +log(t) + exp(-x*x);

    return sign*y;
}
