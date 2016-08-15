#include <iostream>

#include "Algorithm.h"
#include "LikelihoodFreeMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

using namespace MCMC;

int main()
{
    // RNG
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    
    
    // MCMC options
    Options opts;
    gsl_rng_set( r, opts.rng_seed() );
    opts.set_mcmc_trials(1000);
    opts.set_parameter_proposal_sigma(0.1);
    opts.set_observation_noise_sigma(0.1);
    
    double steps = 3;
    using algo_type = Algorithm<LikelihoodFreeMCMC, LangevinDynamics>;
    
    // Sigma
    double log_real_sigma;
    double log_real_sigma_start = -7;
    double log_real_sigma_end = 0;
    double d_step_size = ( log_real_sigma_end - log_real_sigma_start ) / steps;
    std::string subfolder;
    subfolder = "log_real_sigma";
    opts.set_output_subfolder(subfolder);
    for(size_t i=0; i<steps; ++i)
    {
        log_real_sigma = log_real_sigma_start + i*d_step_size;
        std::cout<< "log_real_sigma: " << log_real_sigma << std::endl;
        opts.set_log_real_sigma(log_real_sigma);
        auto algo = algo_type( opts );
        algo.run(r);
    }

}
