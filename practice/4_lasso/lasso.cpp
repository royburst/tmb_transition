#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X)

  // Parameters
 //   PARAMETER(lambda);
  PARAMETER(log_sd);
  PARAMETER_VECTOR(betas);

  // Objective funcction
  Type sd = exp(log_sd);
  Type nll = 0;
  int n = y.size();


  float absbeta = 0.0;
  for(int i = 0; i < betas.size(); i++) absbeta += betas[i] * betas[i];
  absbeta = sqrt(absbeta);

  // priors
  nll = nll - 0.2 * v.sum(absbeta)        // 0.2 is lambda. priors for betas
  nll = nll - log(sd^2)                // prior for sigma^2

  // Linear predictor
  vector<Type> linpred_i(n);
  linpred = X*betas;

  // Probability of data conditional on fixed and random effect values
  for(int i=0; i<n; i++){
    nll -= dnorm(y(i), linpred, sd, true);
  }



  return nll;
}
