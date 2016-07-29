#include "LangevinParams.h"

using namespace MCMC;

void LangevinParamImp::propose_parameters(gsl_rng *r)
{
    int idx;
    double variance = opts.observation_noise_variance();

    for( size_t i=0; i<defining_modes; ++i)
    {
        c_star(i) = gsl_ran_gaussian(r, variance ) + c(i);
    }

    for( size_t i=1; i<cutoff+1; ++i )
        _potential_star.set_mode( i, 0, c_star(i-1) );

    for( size_t i=-cutoff; i<cutoff+1; ++i )
        for( size_t j=1; j<cutoff+1; ++j )
        {
            idx = (cutoff-1) + (i + cutoff)*(2*cutoff + 1) + j;
            _potential_star.set_mode( i, j, c_star( idx ) );
        }
}

double LangevinParamImp::log_transition_density_ratio()
{
    // Symmetric
    return 0;
}

double LangevinParamImp::log_prior_ratio()
{
    double log_total = 0;
    double var = opts.parameter_proposal_variance();
    double normal_constant = 0.5 / var;

    for( size_t i=0; i<defining_modes; ++i)
        log_total -= normal_constant*( pow( std::abs(c_star(i)), 2) -  pow( std::abs(c(i)), 2) );

    return log_total;
}

void LangevinParamImp::accept()
{
    _potential = _potential_star;  
    c_star = c;
}

void LangevinParamImp::store_chain(int n)
{
    for( size_t i=0; i<defining_modes; ++i)
    {
        constants_chain(n,i) = c(i);
    }
}
