#include "LangevinDynamics.h"
#include <iostream>
#include <fstream>

ComplexType LangevinDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double variance = opts.parameter_proposal_variance();
    return ComplexType(gsl_ran_gaussian(r, variance ), gsl_ran_gaussian(r, variance)) + c;
}

double LangevinDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    // Symmetric
    return 0;
}

double LangevinDynamics::log_prior(ParameterType &c)
{
    double log_total = 0;
    double var = opts.parameter_proposal_variance();
    double normal_constant = 0.5 / var;

    int c_dim = parameter_dimension();

    for( size_t i=0; i<c_dim; ++i)
        log_total -= normal_constant*( pow( std::abs(c(i)), 2) );

    return log_total;
}

double LangevinDynamics::log_path_likelihood(  PathType &x, 
                                               ParameterType &c, 
                                               CoarsePathType &y,
                                               int k, int l, int m )
{
    double log_total = 0;
    FourierSeries V(_cutoff);
    V.set_modes(c);
    Vector2d xt;
    Vector2d xtminus1;
    Vector2d gradvtminus1;

    if( m == 0 )
    {
        log_total -= _o_const * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= _o_const * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    {
        xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
        gradvtminus1 = V.grad( xtminus1 );
        log_total -= _d_const * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*_dt, 2);
        log_total -= _d_const * pow( xt(1)-xtminus1(1)+gradvtminus1(1)*_dt, 2);
    }
    
    return log_total;
}

void LangevinDynamics::forward_sim( gsl_rng *r, 
                                    ParameterType &c,
                                    int k, int l, int m, 
                                    PathType &out )
{
    int    M = opts.extra_data_ratio();
    double random_noise_x;
    double random_noise_y;

    FourierSeries V(_cutoff);
    V.set_modes(c);
    
    random_noise_x = gsl_ran_gaussian(r, _d_variance);
    random_noise_y = gsl_ran_gaussian(r, _d_variance);  

    Vector2d xtminus1;
    Vector2d gradtminus1;
    
    if( m != 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);

    gradtminus1 = V.grad( xtminus1 );

    out(k, l, m, 0 ) = xtminus1(0)-gradtminus1(0)*_dt + random_noise_x;
    out(k, l, m, 1 ) = xtminus1(1)-gradtminus1(1)*_dt + random_noise_y;

}

void LangevinDynamics::output_file_timeseries(ParameterChainType &ccc)
{
    std::cout<<"LangevinDynamics::output_file_timeseries"<<std::endl;
    size_t N  = opts.mcmc_trials();

    std::ofstream mcmc_file;
    mcmc_file.open ("output/langevin_mcmc_timeseries.txt");

    int defining_modes = parameter_dimension(opts);

    for( size_t n=0; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(ccc(n,m)) << "\t" << std::imag(ccc(n,m)) << "\t";

        mcmc_file << std::endl;
    }
    mcmc_file.close();
} 
