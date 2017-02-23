#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(n_data);
  DATA_INTEGER(n_factors);
  DATA_IVECTOR(j);
  DATA_VECTOR(y);

  // Parameters
  PARAMETER(X0);
  PARAMETER(log_SD0);
  PARAMETER(log_SDJ);
  PARAMETER_VECTOR(J);

  // Objective funcction
  Type jnll = 0;

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm( y(i), X0 + J(Factor(i)), exp(log_SD0), true );
  }

  // Probability of random coefficients
  for( int i=0; i<n_factors; i++){
    jnll -= dnorm( J(i), Type(0.0), exp(log_SDJ), true );
  }

  // Reporting
  Type SDJ = exp(log_SDJ);
  Type SD0 = exp(log_SD0);

  REPORT(SDJ);
  REPORT(SD0);
  REPORT(J);
  REPORT(X0) ;

  ADREPORT(SDJ);
  ADREPORT(SD0);
  ADREPORT(J);
  ADREPORT(X0);

  return jnll;
}
