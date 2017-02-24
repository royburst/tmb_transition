#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);

  // Parameters
  PARAMETER(log_sd);
  PARAMETER(rho);

  int n_data = y.size();

  // Objective funcction
  Type nll = 0.0;

  //nll -= dnorm(y(0), 0.0, exp(log_sd), true);
  for(int i=1; i<n_data; i++){
    nll -= dnorm(y(i), rho*y[i-1], exp(log_sd)^(1/2), true);
  }

  // Reporting
  Type sd = exp(log_sd);

  REPORT(sd);
  ADREPORT(sd);

  return nll;
}
