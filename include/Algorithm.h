#ifndef ALGORITHM_H
#define ALGORITHM_H

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

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
        Algorithm( Options<Dynamics>& o ) : _opts(o)
        {
            algo_scheme = AlgoType<Dynamics>(o);
        }
        ~Algorithm(){};
        void run( gsl_rng *r );

        Options<Dynamics> opts(){ return _opts; }

    private:
        Options<Dynamics> _opts;
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

#endif
