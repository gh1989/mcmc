#include <fstream>
#include <iostream>

#include "Algorithm.h"
#include "ParticleMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

/* No diffusion constant inference:
./bin/LFLangevin -K 25 -N 10000 -o 0.0002 -p 0.01 -c 0.001 -l -6.8
*/

using namespace MCMC;

int main( int argc, char *argv[] )
{
    Options opts( argc, argv );

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
   
    // Testing SMC
    ParticleMCMC< LangevinDynamics > pmcmc( opts );
    pmcmc.generate_observations(r);
    std::cout << "Observations generated" << std::endl;
        
    pmcmc.generate_true_trajectory(r);
    std::cout<< "True trajectory generated..." << std::endl;
    
    pmcmc.smc(r);
    std::cout<< "Ran SMC algo." << std::endl;
    
    
    
}
