#include <fstream>
#include <iostream>

#include "Algorithm.h"
#include "LikelihoodFreeMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

/*

Drift:

Good.
./bin/LFLangevin -N 5000 -Q 10 -K 50 -M 1 -P 32 -o 0.01 -p 0.0001 -R 123 -c 0.1 -i 0 -B 0

Better.
./bin/LFLangevin -N 5000 -Q 10 -K 50 -M 1 -P 32 -o 0.01 p 0.0001 -R 123 -c 0.04 -i 0 -B 0

Diffusion:
./bin/LFLangevin -N 5000 -K 1 -M 1 -P 8 -o 0.1 -p 0.01 -d 10.0 -i 1 -D 0 -B 1000

K tuning (k=20 is good):
./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.001 -R 123 -c 0.05 -i 0 -B 0

P tuning
./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.001 -c 0.05 -i 0 -B 0

M tuning
./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.001 -c 0.05 -i 0 -B 0

o tuning
./bin/LFLangevin -N 100000 -K 20 -M 1 -P 10 -o 0.01 p 0.001 -c 0.05 -i 0 -B 0



*/

using namespace MCMC;

int main( int argc, char *argv[] )
{
    Options opts( argc, argv );

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
   
    Algorithm<LikelihoodFreeMCMC, LangevinDynamics> algo( opts );
    
    algo.run(r);
    
}
