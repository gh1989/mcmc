#include "BridgeDynamics.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

ComplexType BridgeDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    double random_real_part = gsl_ran_gaussian(r, proposal_sigma);
    double random_imag_part = gsl_ran_gaussian(r, proposal_sigma);
    return ComplexType(random_real_part, random_imag_part) + c;
}

double BridgeDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    return 0;
}

double BridgeDynamics::log_prior(ParameterType &c)
{
    // Changed for one parameter model.
    double log_total = 0;
    double parameter_sigma = _opts.parameter_proposal_sigma();
    double exponential_constant = 0.5 / (parameter_sigma * parameter_sigma);
    int c_dim = parameter_dimension();

    //for( size_t i=0; i<c_dim; ++i)
    log_total -= exponential_constant*( pow( std::abs(c(c_dim-1)), 2) );

    return log_total;
}

// Poorly named: this is just the probability of p( x_{l*M*dt+m*dt} | x_{l*M*dt+m*dt - dt}, c, log_sigma, y )
// unless m==0 then its just the conditional prob p( x_{l*M*dt} | y_{l*M*dt}, c )
double BridgeDynamics::log_path_likelihood(  PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m )
{
    double log_total = 0;
    double dt = _opts.trajectory_path_delta();
    double obs_sigma = _opts.observation_noise_sigma();
    FourierSeries V(_cutoff);
    V.set_modes( c );

    Vector2d xt;
    Vector2d xtminus1;
    Vector2d gradvtminus1;
    
    // dX_t = -grad( V ) dt + sqrt(2D) dW_t
    // D = Exp( log_sigma )
    
    double exponential_constant_path = 0.25 / ( exp(2*log_sigma)*dt );
    double exponential_constant_obs = 0.5 / ( obs_sigma*obs_sigma );
    
    if( m == 0 )
    {        
        // This is always 0... since set starting points without noise currently.
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    {
        xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
        gradvtminus1 = V.grad( xtminus1 );
        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2);
        log_total -= exponential_constant_path * pow( xt(1)-xtminus1(1)+gradvtminus1(1)*dt, 2);
    }
  
    return log_total;
}

void BridgeDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out )
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);
    
    FourierSeries V(_cutoff);
    V.set_modes(c);
    
    double bridge_sigma;
    double delta_m = 1.0 - m*dt;
    double o = _opts.observation_noise_sigma();
    double s = exp( log_sigma );
    double sigma2tilde = pow( s, 2.0 ) / ( pow( s, 2.0 )*delta_m + o );
    
    Vector2d xtminus1;
    Vector2d drift_vector;
    Vector2d end_point;
    Vector2d gradtminus1;
    
    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);

    end_point << y(k,l+1,0), y(k,l+1,1);
    gradtminus1 = V.grad( xtminus1 );
    
    drift_vector = -gradtminus1 + sigma2tilde*( end_point - xtminus1 + gradtminus1 * delta_m );
    bridge_sigma = sqrt( pow( s, 2.0 ) * (1 - sigma2tilde*dt ) );
 
    out(k, l, m, 0 ) = xtminus1(0) + drift_vector(0)*dt + root_2dt*gsl_ran_gaussian(r, bridge_sigma);
    out(k, l, m, 1 ) = xtminus1(1) + drift_vector(1)*dt + root_2dt*gsl_ran_gaussian(r, bridge_sigma);  
}

// Poor name, this is the probability p^( x^k_{l*M*dt+m*dt} | x^k_{l*M*dt+m*dt - dt}, c, log_sigma, y )
double BridgeDynamics::log_bridge_likelihood( PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m)
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);
    
    FourierSeries V(_cutoff);
    V.set_modes(c);
    
    double bridge_variance;
    double delta_m = 1.0 - m*dt;
    double o = _opts.observation_noise_sigma();
    double s = exp( log_sigma );
    double sigma2tilde = pow( s, 2.0 ) / ( pow( s, 2.0 )*delta_m + o );
    
    Vector2d xtminus1;
    Vector2d drift_vector;
    Vector2d end_point;
    Vector2d gradtminus1;
    
    // It must be that m > 0:
    xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
    end_point << y(k,l+1,0), y(k,l+1,1);
    gradtminus1 = V.grad( xtminus1 );
    
    drift_vector = -gradtminus1 + sigma2tilde*( end_point - xtminus1 + gradtminus1 * delta_m );
    bridge_variance = pow( s, 2.0 ) * (1 - sigma2tilde*dt );
 
    double log_total = 0;
    double constant = 0.25 / bridge_variance / dt;
    log_total -= constant*pow( x(k,l,m,0) - xtminus1(0) - drift_vector(0)*dt, 2.0 );
    log_total -= constant*pow( x(k,l,m,1) - xtminus1(1) - drift_vector(1)*dt, 2.0 );
      
    return log_total;
    
}

void BridgeDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file)
{
    int defining_modes = parameter_dimension(_opts);
    
    // headers
    mcmc_file << "#";
    for( size_t m=0; m<defining_modes; ++m)             
        mcmc_file << m << "\t" << m << "\t";
    mcmc_file << "sigma" << std::endl;
    
    // data
    size_t N = _opts.mcmc_trials();
    size_t burn_in = _opts.burn();
    for( size_t n=burn_in; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(ccc(n,m)) << "\t" << std::imag(ccc(n,m)) << "\t";
        mcmc_file << sss(n) << std::endl;
    }
    mcmc_file.close();
} 
