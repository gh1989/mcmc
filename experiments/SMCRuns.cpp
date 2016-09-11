#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include "ParticleMCMC.h"
#include "LangevinDynamics.h"
#include "Options.h"
#include "Particle.h"

using namespace std;

int main( int argc, char *argv[] )
{
    Options opts( argc, argv );

    size_t N = opts.mcmc_trials();
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, opts.rng_seed() );
   
    // Testing SMC
    ParticleMCMC< LangevinDynamics > pmcmc( opts );
    pmcmc.generate_observations(r);
    pmcmc.generate_true_trajectory(r);   
    pmcmc.propose(r);
    double log_marginal;
    
    std::ofstream endpoints_file;
    
    unsigned long int milliseconds_since_epoch = 
    std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    
    std::string filename = "output/SMC/";
    filename += opts.output_subfolder();
    filename += "/SMCEndPoints_";
    filename += std::to_string( milliseconds_since_epoch );
    filename += ".txt";
    endpoints_file.open(filename);
    
    endpoints_file << "#P " << opts.path_length() << std::endl;
    endpoints_file << "#M " << opts.extra_data_ratio() << std::endl;
    endpoints_file << "#Q " << opts.number_particles() << std::endl;
    
    for( size_t n=0; n<N; ++n )
    {
        std::cout << "Running SMC..." << n << std::endl;
        pmcmc.propose(r);
        pmcmc.smc(r);
        pmcmc.log_particle_end_points(endpoints_file);
    }    
    
    endpoints_file.close();
}
