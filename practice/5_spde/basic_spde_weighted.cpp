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
  // Indices
  DATA_INTEGER( n_i ); // number of datapts in space time
  DATA_INTEGER( n_x ); // number of mesh pts in space mesh
  DATA_INTEGER( n_t ); // number of time periods
  DATA_INTEGER( n_p ); // number of covariate columns in design matrix X

  // Data
  //  DATA_IVECTOR( x_s );  // <DEPRECATED> association of each cluster with a given vertex
  DATA_VECTOR( c_i );   // obs deaths per cluster
  DATA_VECTOR( Exp_i ); // exposure mths per cluster
  DATA_IVECTOR( s_i );  // observation specific id
  DATA_IVECTOR( t_i );  // # time periods
  DATA_IVECTOR( w_i );  // weights for observations
  DATA_MATRIX( X_xp );  // covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(G0); // used to make gmrf precision
  DATA_SPARSE_MATRIX(G1); // used to make gmrf precision
  DATA_SPARSE_MATRIX(G2); // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Options
  DATA_VECTOR( options ) // basically for if/else statments on the cpp side for model flexibility

  // Fixed effects
  PARAMETER_VECTOR(alpha); // fixed effect coefs
  PARAMETER(log_tau_E);    // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(log_kappa);    // log of INLA kappa - related to spatial correlation and range
  PARAMETER(rho);          // temporal autocorrelation in AR1

  // Random effects
  PARAMETER_ARRAY(epsilon); // GMRF Sp-Tm random effects
  //  PARAMETER_ARRAY(sp); // <DEPRECATED>


  // objective function -- joint negative log-likelihood
  using namespace density;
  vector<Type> jnll_comp(3);
  jnll_comp[0] = Type(0);
  jnll_comp[1] = Type(0);
  jnll_comp[2] = Type(0);
  //  parallel_accumulator<Type> jnll_comp[1](this);
  max_parallel_regions = omp_get_max_threads();
  printf("This is thread %d\n", max_parallel_regions);

  // Spatial parameters
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;
  Type rho_trans = log((1+rho)/(1-rho));


  // Objects for derived values
  vector<Type> linear_x(n_x);          // main effect X*alpha
  matrix<Type> Epsilon_xt(n_x, n_t);   // matrix of epsilons scaled by tau
  vector<Type> epsilon_vec(n_x * n_t); // Epsilon_xt unlisted into a vector for easier matrix multiplication
  vector<Type> proj_epsilon(n_i);      // value of gmrf at data points

  ///////////////////////////////////
  // EVALUATE THE JOINT LIKELIHOOD //
  ///////////////////////////////////

  // Priors
  if(options[0] == 1) {
   jnll_comp[2] -= dnorm(log_tau_E, Type(0.0), Type(1.0), true);  // N(0,1) prior for log_tau
   jnll_comp[2] -= dnorm(log_kappa, Type(0.0), Type(1.0), true);  // N(0,1) prior for log_kappa
   jnll_comp[2] -= dnorm(rho_trans, Type(0.0), Type(2.582), true); // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   for( int i=0; i < alpha.size(); i++){
     // printf("This is alpha %d\n", i);
     jnll_comp[2] -= dnorm(alpha(i), Type(0.0), Type(100), true); // N(0, sqrt(1/.0001)) prior for fixed effects.
   }
  }

  // Probability of Gaussian-Markov random fields (GMRFs)
  if(n_t > 1 ){ // if more than 1 time period, then AR1 in time
    PARALLEL_REGION jnll_comp[0] += SCALE(SEPARABLE(AR1(rho),GMRF(Q)),1/exp(log_tau_E))(epsilon);
  } else { // otherwise just spatial
    PARALLEL_REGION jnll_comp[0] += SCALE(GMRF(Q),1/exp(log_tau_E))(epsilon);
  }


  // Transform GMRFs and make vector form
  for(int x=0; x<n_x; x++){
    for(int t=0; t<n_t; t++){
      Epsilon_xt(x,t) = epsilon(x,t) / exp(log_tau_E); // scale it // we probably don't need this anymore
      epsilon_vec[(x + n_x*t)] = Epsilon_xt(x,t); // put in a vector
    }
  }


  // Project from mesh points to data points
  proj_epsilon = Aproj * epsilon_vec.matrix();

  // Likelihood contribution from observations
  linear_x = X_xp * alpha.matrix();
  vector<Type> mrprob(n_i);
  for (int i=0; i<n_i; i++){
    // THIS LINE CAUSED THE SCALING BUG:    Epsilon_xt(x_s(s_i(i)),t_i(i)) = epsilon(x_s(s_i(i)),t_i(i))/exp(log_tau_E);
    // mrprob(i) = linear_x(x_s(s_i(i))) + Epsilon_xt(x_s(s_i(i)),t_i(i)); // <DEPRECATED> to accont for datapts not at mesh locs
    mrprob(i) = linear_x(i) + proj_epsilon(i); // linear part plus projected sp-tm gp
    if( !isNA(c_i(i)) ){
       PARALLEL_REGION jnll_comp[1] -= dbinom( c_i(i), Exp_i(i), invlogit(mrprob(i)), true )*w_i(i);
    // TEST   PARALLEL_REGION jnll_comp[1] -= dpois(  c_i(i), invlogit(mrprob(i)) * Exp_i(i), true);
    }
  }
  Type jnll = jnll_comp.sum();


  // Diagnostics
  // REPORT( jnll_comp );
  // REPORT( jnll );
  // REPORT( SigmaE );
//  REPORT( alpha );
  ADREPORT( alpha );
//  REPORT( Epsilon_xt );
  ADREPORT( Epsilon_xt );


  return jnll;
}
