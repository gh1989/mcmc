#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fstream>
#include <iostream>

#include <fstream>
#include <string>

#include "Options.h"
#include "LangevinDynamics.h"

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>


/*
This is a test to work out the distribution of trajectory termination.
*/

using namespace MCMC;

int main(int argc, char *argv[])
{
    auto o  = Options( argc, argv );
    auto ld = LangevinDynamics(o);

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, o.rng_seed() );

    size_t N = o.mcmc_trials();
    size_t K = o.parallel_paths(); 
    size_t M = o.extra_data_ratio();
    size_t L = o.path_length();

    LangevinDynamics::PathType paths(K,L,M,2);
    LangevinDynamics::CoarsePathType observed(K,L,2);

    for(size_t k=0; k<K; ++k)
    {
        observed(0, 0, 0 ) = 0.0;
        observed(0, 0, 1 ) = 0.0;
    }

    auto c = LangevinDynamics::default_parameters(o);
    double log_sigma = log( 0.001 );

    Tensor<double, 3> end_points(N,K,2);

    std::string ts_file_name = "output/end_points.txt";
    std::ofstream ts_file;
    ts_file.open(ts_file_name);

    for(size_t n=0; n<N; ++n)
    {
    for(size_t k=0; k<K; ++k)
    {
        for(size_t l=0; l<L; ++l)
        {
            for(size_t m=0; m<M; ++m)
            {
                ld.forward_sim(r, c, log_sigma, observed, 0, l, m, paths);   
            }
        }
        end_points(n, k, 0 ) = paths(k,L-1,0,0);
        end_points(n, k, 1 ) = paths(k,L-1,0,1);
    }    
    }

    for(size_t n=0; n<N; ++n )
    {
        for(size_t k=0; k<K; ++k)
        {
            ts_file << end_points(n,k,0) << "\t" << end_points(n,k,1) << std::endl;
        }    
    }   

    ts_file.close();

}
