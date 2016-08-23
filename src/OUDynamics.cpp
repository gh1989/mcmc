#include "OUDynamics.h"
#include <iostream>
#include <fstream>

using namespace MCMC;


double OUDynamics::sample_transition_density(gsl_rng *r, double c )
{
    double sigma = _opts.parameter_proposal_sigma();
    return gsl_ran_gaussian(r, sigma ) + c;
}

double OUDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    // Symmetric
    return 0;
}

double OUDynamics::log_prior(ParameterType &c)
{
    double sigma = _opts.parameter_proposal_sigma();
    double exponential_constant = 0.5 / sigma / sigma;
    return -exponential_constant*( pow( c(0), 2) );
}

double OUDynamics::log_path_likelihood( PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m )
{
    double log_total = 0;
    double dt = _opts.trajectory_path_delta();
    double obs_sigma = _opts.observation_noise_sigma();
    double exponential_constant_path = 0.25/(exp(log_sigma)*exp(log_sigma)*dt);
    double exponential_constant_obs  = 0.5/(obs_sigma*obs_sigma);
    
    Vector2d xt;
    Vector2d xtminus1;

    if( m == 0 )
    {
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    { 
        xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 ); 

        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(0)+c(0)*xtminus1(0)*dt, 2);
        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(1)+c(1)*xtminus1(1)*dt, 2);        
    }
    return log_total;
}

void OUDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out )
{
    int M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double exponential_constant_path = 0.25/(exp(log_sigma)*exp(log_sigma)*dt);
    double root_2dt = sqrt( 2*dt );

    Vector2d xtminus1;

    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);

    out(k, l, m, 0 ) = xtminus1(0)-c(0)*xtminus1(0)*dt + root_2dt*exp( 0.5*log_sigma )*gsl_ran_gaussian(r, 1 );
    out(k, l, m, 1 ) = xtminus1(1)-c(0)*xtminus1(1)*dt + root_2dt*exp( 0.5*log_sigma )*gsl_ran_gaussian(r, 1 );

}


void OUDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file)
{

    std::cout<<"OUDynamics::output_file_timeseries"<<std::endl;
    size_t N  = _opts.mcmc_trials();
    size_t burn_in = _opts.burn();
    for( size_t n=burn_in; n<N; ++n )
    {
        mcmc_file << ccc(n,0);
        mcmc_file << "\t";
        mcmc_file << sss(n);
        mcmc_file << std::endl;
    }

}
