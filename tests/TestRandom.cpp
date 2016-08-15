#include <iostream>
#include <fstream>
#include <string>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

/*
The purpose of these tests is to ensure the distributions of the random number
generators is as expected. 
*/

int main()
{
    int seed = 123;
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, seed );
    
    string filename = "output/test_gsl_ran_gaussian.txt";
    ofstream mcmc_file;
    mcmc_file.open(filename);
    
    // First test the gsl_ran_gaussian function.
    double sigma = 0.25;
    for( size_t i=0; i<10000; ++i )
    {
        mcmc_file << gsl_ran_gaussian(r, sigma) << std::endl;
    }
    
    mcmc_file.close();
    
}