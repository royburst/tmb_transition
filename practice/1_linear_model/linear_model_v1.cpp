// Just always declare your functions like they are in the 4 lines below
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
  Type jnll = 0;            // starting point for nll (the 'objective function')
                            // Use Type instead of float or double
  int n_data = y_i.size();  // data length value to loop over

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm(y_i(i), mean, sd, true); // subtract on the value t
  }

  // I like loops and i think they dont make things slower. This vector approach
  // would work as well, I think: jnll = -sum(dnorm(y_i,mu,sigma,true));

  // Reporting
  return jnll;
}
