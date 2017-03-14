#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
//  parallel_accumulator<Type> nll(this);
  using namespace density;

  Type nll = 0;
  max_parallel_regions = omp_get_max_threads();

  for(int i=0;i<x.size();i++){
    PARALLEL_REGION nll -= dnorm(Y[i],a+b*x[i],exp(logSigma),true);
  }

  return nll;
}
