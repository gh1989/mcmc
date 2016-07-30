#include "Algorithm.h"
#include <iostream>

using namespace MCMC;

void Algorithm::run( gsl_rng *r )
{
    double log_a;
    double log_u;
    double acceptance_rate = 0.0;    

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
        {
           acceptance_rate += 1.0;         
           std::cout<<"Accept ("<<double(acceptance_rate)/double(n)<<")"<<std::endl;          
           algo_scheme.accept();
        }
        else
        {
            std::cout<<"Reject ("<<double(n-acceptance_rate)/double(n)<<")"<<std::endl;     
        }
        algo_scheme.store_chain(n);
        
    }

    std::cout<< "Accepted: " << double(acceptance_rate)/double(N) << std::endl;
    algo_scheme.finish();
}
