#include "Algorithm.h"
#include <iostream>

using namespace MCMC;

void Algorithm::run( gsl_rng *r )
{
    double log_a;
    double log_u;
    
    int N = _opts.mcmc_trials();

    std::cout<<"Starting MCMC algorithm. N = "<< N << std::endl;

    algo_scheme.generate_true_trajectory(r);
    algo_scheme.generate_observations(r);
    
    for(size_t n=0; n<N; ++n)
    {   
        algo_scheme.propose_parameters(r);
        algo_scheme.propose_trajectory(r);
        
        log_u = log( gsl_rng_uniform(r) );
        log_a = algo_scheme.log_acceptance_probability();   

        if( log_u < log_a )
            algo_scheme.accept();
        
        algo_scheme.store_chain(n);

    }
    algo_scheme.finish();
}
