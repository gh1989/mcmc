#include "OUParams.h"

void OUParamImp::propose_parameters(gsl_rng *r)
{
    double variance = opts.parameter_proposal_variance();
    _c_star = gsl_ran_gaussian(r, variance) + _c;   
}

double OUParamImp::log_transition_density_ratio()
{   
    // Symmetric.
    return 0;
}

double OUParamImp::log_prior_ratio()
{
    double variance = opts.parameter_proposal_variance();
    return -(0.5 / variance ) * ( pow( _c_star, 2 ) - pow( _c, 2) ); 
}
                      

void OUParamImp::accept()
{
    _c = _c_star;
    sigma = sigma_star;
}

void OUParamImp::store_chain(int n)
{
    constants_chain(n) = _c;
    diffusion_chain(n) = sigma;
}

void OUParamImp::output_file_timeseries()
{
    size_t N  = opts.mcmc_trials();

    std::ofstream ou_mcmc_file;
    ou_mcmc_file.open ("output/ou_mcmc_timeseries.txt");

    for( size_t n=0; n<N; ++n )
    {
        ou_mcmc_file << constants_chain(n) << std::endl;
    }
    ou_mcmc_file.close();
}

