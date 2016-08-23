#include <fstream>
#include <iostream>

#include "Algorithm.h"
#include "Options.h"
#include "LangevinDynamics.h"

/*
 *
 * An experiment to use the LF-MCMC algorithm to infer a diffusion coefficient
 * for the SDE $dX_t = \sigma dW_t$.
 *
 */

using namespace MCMC;

int main( int argc, char *argv[] )
{
    // Collect options from command line arguments.
    Options opts( argc, argv );
  
    // The random number generator.
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
   
    // Infer the diffusion coefficient.
    opts.set_infer_diffusion_parameters(true);

    // We will not infer drift parameters here. 
    opts.set_infer_drift_parameters(false);

    // Instantiate an LF-MCMC with Langevin dynamics algo then run it.
    Algorithm<LikelihoodFreeMCMC, LangevinDynamics> algo( opts );
    algo.run(r);
 
}
