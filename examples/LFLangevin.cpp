#include <iostream>

#include "Algorithm.h"
#include "LikelihoodFreeMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

/*

OUDynamics hard-coded constant: true_c = 1.70;

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
    Options<LangevinDynamics> opts( argc, argv );

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );

    LangevinDynamics::ParameterType real_parameters( 4 );
    real_parameters(3) = ComplexType( 0.5, -0.5 );
    opts.set_parameters( real_parameters ); 

    Algorithm<LikelihoodFreeMCMC, LangevinDynamics> algo( opts );
    
    algo.run(r);
    
}
