// Space time
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Indices
  DATA_INTEGER( n_i );
  DATA_INTEGER( n_x );
  DATA_INTEGER( n_t );
  DATA_INTEGER( n_p );

  // Data
  DATA_IVECTOR( x_s );
  DATA_VECTOR( c_i );
  DATA_VECTOR( Exp_i );
  DATA_IVECTOR( s_i );
  DATA_IVECTOR( t_i );
  DATA_MATRIX( X_xp );

  // SPDE objects
  DATA_STRUCT(spde,spde_t);

  // Fixed effects
  PARAMETER_VECTOR(alpha);
  PARAMETER(log_tau_E);
  PARAMETER(log_kappa);

  // Random effects
  PARAMETER_ARRAY(epsilon);
//  PARAMETER_ARRAY(sp);


  // objective function -- joint negative log-likelihood

  Type jnll = 0;
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
//    parallel_accumulator<Type> jnll(this);
  max_parallel_regions = omp_get_max_threads();
  printf("This is thread %d\n", max_parallel_regions);

  // Spatial parameters
  Type kappa = exp(log_kappa);
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
//  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;
//  Type rho_trans = log((1+rho)/(1-rho));
  SparseMatrix<Type> Q = Q_spde(spde,kappa);

  // Objects for derived values
  vector<Type> linear_x(n_x);
  matrix<Type> Epsilon_xt(n_x, n_t);

  // Priors
  // jnll_comp(3) -= dnorm(log_tau_E, Type(0.0), Type(1.0), true);  // N(0,1) prior for log_tau
  // jnll_comp(3) -= dnorm(log_kappa, Type(0.0), Type(1.0), true);  // N(0,1) prior for log_kappa
  // jnll_comp(3) -= dnorm(rho_trans, Type(0.0), Type(2.582), true); // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
  // jnll_comp(3) -= dnorm(alpha, Type(0.0), Type(100), true); // N(0, sqrt(1/.0001)) prior for fixed effects.
  // NOTE: in INLA the intercept is given a prior N(mean=0, prec=0) which I don't know how to code here


  // Probability of Gaussian-Markov random fields (GMRFs)
//   jnll += GMRF(Q)(sp);
   if(n_t > 1 ){
     PARAMETER(rho);
     jnll_comp[0] += SEPARABLE(AR1(rho),GMRF(Q))(epsilon);
   } else {
     jnll_comp[0] += GMRF(Q)(epsilon);
   }


  // Transform GMRFs
  for(int x=0; x<n_x; x++){
    for( int t=0; t<n_t; t++){
      Epsilon_xt(x,t) = epsilon(x,t) / exp(log_tau_E);
    }
  }

  // Likelihood contribution from observations
  linear_x = X_xp * alpha.matrix();
  vector<Type> mrprob(n_i);
  for (int i=0; i<n_i; i++){
    mrprob(i) = linear_x(x_s(s_i(i))) + Epsilon_xt(x_s(s_i(i)),t_i(i));
    if( !isNA(c_i(i)) ){
       PARALLEL_REGION jnll_comp[1] -= dbinom( c_i(i), Exp_i(i), invlogit(mrprob(i)), true );
    }
  }
  jnll = jnll_comp.sum();


  // Diagnostics
  // REPORT( jnll_comp );
  // REPORT( jnll );
  // REPORT( SigmaE );
  REPORT( alpha );
  ADREPORT( alpha );
  REPORT( Epsilon_xt );
  ADREPORT( Epsilon_xt );


  return jnll;
}
