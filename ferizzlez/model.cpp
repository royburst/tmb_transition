// ///////////////////////////////////////////////////
// Roy Burstein and Aaron Osgood-Zimmerman
// August 2017
// Template file for space-time-Z GPR model.
// Used for fitting IHME Geospatial MBG models
// ///////////////////////////////////////////////////

// ///////////////////////////////////////////////////
// NOTES:
// 1. Type `Type` is a special TMB type that should be used for all variables, except `int` can be used as well
// 2. In our nomenclature, Z is a third interaction (ie age) which defaults to AR1
// 3. Requires same space mesh for all time-Z points
// 4. Anything in the density namespace (ie SEPARABLE) returns the negative log likelihood and is thus added to accumulator
//    also, other density function such as dnorm and dbinom return positive log likelihood and are thus subtracted away.
// 5. ref https://github.com/nmmarquez/re_simulations/blob/master/inla/sta.cpp
//        https://github.com/nmmarquez/re_simulations/blob/master/inla/SPDEAR1AR1.R
// ///////////////////////////////////////////////////

// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;


// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// objective function (ie the likelihood function for the model), returns the evaluated negative log likelihood
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ////////////////////////////////////////////////////////////////////////////
  // INPUTS
  // ////////////////////////////////////////////////////////////////////////////

  // Indices
  DATA_INTEGER(num_i);       // number of datapts in space-time-Z (aka STZ)
  DATA_INTEGER(num_s);       // number of mesh pts in space mesh
  DATA_INTEGER(num_t);       // number of time periods
  DATA_INTEGER(num_z);       // number of Z groups

  // Data (each, excpect for X_ij is a vector of length num_i)
  DATA_VECTOR(y_i);          // obs successes per binomial experiment at point i (aka cluster)
  DATA_VECTOR(n_i);          // trials per cluster
  DATA_IVECTOR(t_i);         // time period of the data point
  DATA_IVECTOR(w_i);         // weights for observations
  DATA_MATRIX(X_ij);         // covariate design matrix (num_i by number of fixed effects matrix)

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations

  // Options
  DATA_VECTOR(options)       // boolean vector of options to be used to select different models/modelling options:
                             // 0: Include priors. All are default settings right now
                             // 1:
                             // 2:

  // Parameters
  PARAMETER_VECTOR(alpha_j); // fixed effect coefs, including intercept as first index
  PARAMETER(logtau);         // log of INLA tau param (precision of space-time covariance mat)
  PARAMETER(logkappa);       // log of INLA kappa - related to spatial correlation and range
  PARAMETER(trho);           // temporal autocorrelation parameter for AR1, natural scale
  PARAMETER(zrho);           // Z autocorrelation parameter for AR1, natural scale

  // Random effects
  PARAMETER_ARRAY(Epsilon_stz); // Random effects for each STZ mesh location. Should be 3D array of dimensions num_s by num_t by num_z

  printf("Epsilon_stz dimensions: %d\n", Epsilon_stz.size());

  // ////////////////////////////////////////////////////////////////////////////
  // LIKELIHOOD
  // ////////////////////////////////////////////////////////////////////////////

  // Define the joint-negative log-likelihood as a parallel_accumulator
  // this allows us to add or subtract numbers to the object in parallel
  parallel_accumulator<Type> jnll(this);


  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(logkappa, logtau, M0, M1, M2);


  // Make transformations of some of our parameters
  Type range     = sqrt(8.0) / exp(logkappa);
  Type sigma     = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) * exp(2.0 * logkappa));
  Type trho_trans = log((1.0 + trho) / (1.0 - trho));
  Type zrho_trans = log((1.0 + zrho) / (1.0 - zrho));


  // Define objects for derived values
  vector<Type> fe_i(num_s);                         // main effect X_ij %*% t(alpha_j)
  vector<Type> epsilon_stz(num_s * num_t * num_z);  // Epsilon_stz unlisted into a vector for easier matrix multiplication
  vector<Type> projepsilon_i(num_i);                // value of gmrf at data points
  vector<Type> prob_i(num_i);                       // Logit estimated prob for each point i


  // Prior contribution to likelihood. Values are defaulted (for now). Only run if options[0]==1
  if(options[0] == 1) {
   jnll -= dnorm(logtau,    Type(0.0), Type(1.0),   true);  // N(0,1) prior for logtau
   jnll -= dnorm(logkappa,  Type(0.0), Type(1.0),   true);  // N(0,1) prior for logkappa
   if(num_t > 1) {
     jnll -= dnorm(trho_trans, Type(0.0), Type(2.582), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   if(num_z > 1) {
     jnll -= dnorm(zrho_trans, Type(0.0), Type(2.582), true);  // N(0, sqrt(1/.15) prior on log((1+rho)/(1-rho))
   }
   for( int j = 0; j < alpha_j.size(); j++){
     jnll -= dnorm(alpha_j(j), Type(0.0), Type(100), true); // N(0, sqrt(1/.001)) prior for fixed effects.
   }
  }

  // Latent field/Random effect contribution to likelihood.
  // Possibilities of Kronecker include: S, ST, SZ, and STZ
  if (num_t == 1 & num_z == 1)  {
    printf("SPACE  ONLY");
    jnll += GMRF(Q_ss)(epsilon_stz);
  } else if(num_t > 1 & num_z == 1) {
    printf("SPACE-TIME");
    jnll += SEPARABLE(AR1(trho),GMRF(Q_ss))(Epsilon_stz);
  } else if (num_t == 1 & num_z > 1) {
    printf("SPACE-Z");
    jnll += SEPARABLE(AR1(zrho),GMRF(Q_ss))(Epsilon_stz);
  } else if (num_t > 1 & num_z > 1) {
    printf("SPACE-TIME-Z");
    jnll += SEPARABLE(AR1(zrho),SEPARABLE(AR1(trho),GMRF(Q_ss)))(Epsilon_stz);
  }

  // Transform GMRFs and make vector form
  // TODO check indexing
  for(int s = 0; s < num_s; s++){
    for(int t = 0; t < num_t; t++){
      for(int z = 0; z < num_z; z++){
        epsilon_stz[(s + num_s * t + num_t * z)] = Epsilon_stz(s,t,z); // put in a vector
      }
    }
  }

  // Project from mesh points to data points in order to eval likelihood at each data point
  // TODO expand this for Z
  projepsilon_i = Aproj * epsilon_stz.matrix();

  // evaluate fixed effects for alpha_j values
  fe_i = X_ij * alpha_j.matrix();

  // Likelihood contribution from each datapoint i
  for (int i = 0; i < num_i; i++){
    prob_i(i) = fe_i(i) + projepsilon_i(i);
    if(!isNA(y_i(i))){
      jnll -= dbinom( y_i(i), n_i(i), invlogit(prob_i(i)), true ) * w_i(i);
    }
  }

  // Report estimates
  ADREPORT(alpha_j);
  ADREPORT(Epsilon_stz);

  return jnll;
}
