#include "OUDynamics.h"
#include <iostream>
#include <fstream>

using namespace MCMC;

double OUDynamics::sample_transition_density(gsl_rng *r, double c )
{
    double variance = opts.parameter_proposal_variance();
    return gsl_ran_gaussian(r, variance ) + c;
}

double OUDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    // Symmetric
    return 0;
}

double OUDynamics::log_prior(ParameterType &c)
{
    double var = opts.parameter_proposal_variance();
    double normal_constant = 0.5 / var;
    return -normal_constant*( pow( c(0), 2) );
}

double OUDynamics::log_path_likelihood( PathType &x, 
                                        ParameterType &c, 
                                        double sigma,
                                        CoarsePathType &y,
                                        int k, int l, int m )
{
    double log_total = 0;
    double _d_const = 0.5/(sigma*sigma*_dt);
    
    Vector2d xt;
    Vector2d xtminus1;

    if( m == 0 )
    {
        log_total -= _o_const * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= _o_const * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    { 
        xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 ); 

        log_total -= _d_const * pow( xt(0)-xtminus1(0)+c(0)*xtminus1(0)*_dt, 2);
        log_total -= _d_const * pow( xt(0)-xtminus1(1)+c(1)*xtminus1(1)*_dt, 2);        
    }
    return log_total;
}

void OUDynamics::forward_sim( gsl_rng *r, 
                              ParameterType &c,
                              double sigma,
                              int k, int l, int m, 
                              PathType &out )
{
    int M = opts.extra_data_ratio();

    double random_noise_x;
    double random_noise_y;
    double _d_variance = 0.5 / (sigma*sigma*_dt);
    
    random_noise_x = gsl_ran_gaussian(r, _d_variance );
    random_noise_y = gsl_ran_gaussian(r, _d_variance);  

    Vector2d xtminus1;

    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);

    out(k, l, m, 0 ) = xtminus1(0)-c(0)*xtminus1(0)*_dt + random_noise_x;
    out(k, l, m, 1 ) = xtminus1(1)-c(0)*xtminus1(1)*_dt + random_noise_y;

}


void OUDynamics::output_file_timeseries(ParameterChainType &ccc)
{

    std::cout<<"OUDynamics::output_file_timeseries"<<std::endl;
    size_t N  = opts.mcmc_trials();

    std::ofstream mcmc_file;
    mcmc_file.open ("output/ou_mcmc_timeseries.txt");

    for( size_t n=0; n<N; ++n )
    {
        mcmc_file << ccc(n,0);
        mcmc_file << std::endl;
    }
    mcmc_file.close();

}
