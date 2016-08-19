#include <fstream>
#include <iostream>

#include "Algorithm.h"
#include "ParticleMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

/* No diffusion constant inference:
./bin/LFLangevin -K 25 -N 10000 -o 0.0002 -p 0.01 -c 0.001 -l -6.8
./bin/PMCMC -Q 10 -K 10 -p 0.0001 -M 4 -P 4 -o 0.001 -R 123 -i 0
*/

using namespace MCMC;

int main( int argc, char *argv[] )
{
    Options opts( argc, argv );

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
   
    Algorithm<ParticleMCMC, LangevinDynamics> algo( opts );
    
    algo.run(r);
    
}
