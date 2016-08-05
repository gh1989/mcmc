#include <iostream>

#include "Algorithm.h"
#include "SequentialLFMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

using namespace MCMC;

int main( int argc, char *argv[] )
{
    Options<LangevinDynamics> opts( argc, argv );
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
    Algorithm<SequentialLFMCMC, LangevinDynamics> algo( opts );
    algo.run(r);   
}
