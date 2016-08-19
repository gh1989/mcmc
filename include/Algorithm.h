#ifndef ALGORITHM_H
#define ALGORITHM_H

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LikelihoodFreeMCMC.h"
#include "Options.h"

using namespace MCMC;

namespace MCMC
{

template <template <class A> class AlgoType, class Dynamics>
class Algorithm
{
    public:
    
        using algo_type = AlgoType<Dynamics>;
    
        Algorithm( Options& o ) : _opts(o), algo_scheme(o) 
        {
            _opts.print_options(std::cout);
        }
        
        Algorithm(Algorithm&&)=default;
        
        template< template <class X> class U, class V >
        Algorithm& operator=(Algorithm<U, V>&& arg)
        {
        };      
        
        Algorithm(Algorithm&)=default;       
        
        template< template <class X> class U, class V >
        Algorithm& operator=(Algorithm<U, V>& arg)
        {
            swap(arg);
            return *this;
        };     
        
        ~Algorithm() {};
        void run( gsl_rng *r );

        Options opts(){ return _opts; }

    private:
        Options _opts;
        AlgoType<Dynamics> algo_scheme;
};

}

template <template <class A> class AlgoType, class Dynamics>
void Algorithm<AlgoType, Dynamics>::run( gsl_rng *r )
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
        algo_scheme.propose(r);
        log_u = log( gsl_rng_uniform(r) );
        log_a = algo_scheme.log_acceptance_probability();   

        if( log_u < log_a )
        {
           acceptance_rate += 1.0;         
           std::cout<< n << "... Accept ("<<double(acceptance_rate)/double(n+1)<<")"<<std::endl;             
           algo_scheme.accept();
        }
        else
        {
            std::cout<< n <<"... Reject"<<std::endl;           
        }
        
        algo_scheme.store_chain(n);
        
        if (_opts.infer_drift_parameters())
        {
            std::cout << "c" << std::endl;
            std::cout << algo_scheme.current_drift() << std::endl;
        }
        
        if (_opts.infer_diffusion_parameters())
        {
            std::cout << "log_sigma" << std::endl;
            std::cout << algo_scheme.current_log_sigma() << std::endl;
        }
        
    }

    std::cout<< std::endl << "Accepted: " << double(acceptance_rate)/double(N) << std::endl;  
    algo_scheme.finish();
}

#endif
