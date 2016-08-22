#include <iostream>

#include "Algorithm.h"
#include "LikelihoodFreeMCMC.h"
#include "Options.h"
#include "OUDynamics.h"

/*

OUDynamics hard-coded constant: true_c = 1.70;

set _mcmc_trials                    = 1000
set _trajectory_path_delta          = 0.05
set _parallel_paths                 = 500
set _extra_data_ratio               = 1
set _observation_noise_sigma     = 0.01
set _path_length                    = 16
set _parameter_proposal_variance    = 0.05
set _diffusion_coefficient          = 0.001
set _burn                           = 0;
set _rng_seed                       = 0;
set _parameters                     = 1;
*/

using namespace MCMC;

int main( int argc, char *argv[] )
{
    Options opts( argc, argv );

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
 
    Algorithm<LikelihoodFreeMCMC, OUDynamics> algo( opts );
    
    algo.run( r );
    
}


