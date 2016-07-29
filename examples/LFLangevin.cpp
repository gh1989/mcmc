#include <iostream>

#include "Algorithm.h"
#include "LikelihoodFreeScheme.h"
#include "Options.h"
#include "LangevinPathScheme.h"

/*

OUPathScheme hard-coded constant: true_c = 1.70;

set _mcmc_trials                    = 10000
set _trajectory_path_delta          = 0.02
set _parallel_paths                 = 25
set _extra_data_ratio               = 1
set _observation_noise_variance     = 0.1
set _path_length                    = 16
set _parameter_proposal_variance    = 0.1
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
 
    LangevinPathScheme langevin_scheme( opts );
    LikelihoodFreeScheme lf_langevin( opts,  langevin_scheme );
    Algorithm algo( opts, lf_langevin );
    
    algo.run(r);
    
}
