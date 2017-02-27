#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X);

  // Parameters
 //   PARAMETER(lambda);
  PARAMETER(log_sd);
  PARAMETER(L);
  PARAMETER_VECTOR(betas);

  // Objective funcction
  Type sd = exp(log_sd);
  Type nll = 0;
  int n = y.size();

  // get absolute values of the betas
  Type absbeta = 0.0;
  for(int i = 0; i < betas.size(); i++) absbeta += sqrt(pow(betas(i),2));

  // priors
  nll = nll + dunif(L,0.001, 3, true);          // lambda prior
  nll = nll + L * absbeta;             // 0.2 is lambda. priors for betas
  nll = nll + log(pow(sd,2));             // prior for sigma^2

  // Linear predictor
  vector<Type> linpred(n);
  linpred = X*betas;

  // Probability of data conditional on fixed and random effect values
  for(int i=0; i<n; i++){
    nll -= dnorm(y(i), linpred(i), sd, true);
  }



  return nll;
}
