#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y_i);
  DATA_MATRIX(X_ij)

  // Parameters
  PARAMETER(lambda);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(betas);

  // Objective funcction
  Type sd = exp(log_sd);
  Type nll = 0;
  int n_data = y_i.size();


  // Probability of data conditional on fixed and random effect values
  for(int i=0; i<n_data; i++){
    jnll -= dnorm(y(i), X0 + J(j(i)), exp(log_SD0), true);
  }

  // Probability of random coefficients
  for( int i=0; i<n_factors; i++){
    jnll -= dnorm(J(i), Type(0.0), exp(log_SDJ), true);
  }

  // Reporting
  Type SDJ = exp(log_SDJ);
  Type SD0 = exp(log_SD0);


  return jnll;
}
