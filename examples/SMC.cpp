#include <fstream>
#include <iostream>

#include "Algorithm.h"
#include "ParticleMCMC.h"
#include "Options.h"
#include "LangevinDynamics.h"

/* No diffusion constant inference:
./bin/SMC -Q 1000 -K 10 -M 10 -P 4 -o 0.01 -p 0.001 -R 1234 -i 0 -D 1
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
    
    pmcmc.propose(r);
    double log_marginal = pmcmc.smc(r);
    std::cout<< "Ran SMC algo." << std::endl;
    std::cout<< "Log marginal: " << log_marginal << std::endl;
        
}
