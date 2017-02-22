// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y_i);

  // Parameters
  PARAMETER(mean);
  PARAMETER(log_sd);

  // Objective funcction
  Type sd = exp(log_sd);    // ensures that SD is positive
  Type jnll = 0;            // starting point for nll
  int n_data = y_i.size();  // data length value to loop over

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm(y_i(i), mean, sd, true); // subtract on the value t
  }

  // Reporting
  return jnll;
}
