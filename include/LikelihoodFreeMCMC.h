#ifndef LIKELIHOOD_MCMC_H
#define LIKELIHOOD_MCMC_H

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "Dynamics.h"

namespace MCMC
{

template<typename Dynamics_>
class LikelihoodFreeMCMC
{
public:

    typedef typename Dynamics_::ParameterType ParameterType;
    typedef typename Dynamics_::ParameterChainType ParameterChainType;
    typedef typename Dynamics_::PathType PathType;
    typedef typename Dynamics_::CoarsePathType CoarsePathType;

    LikelihoodFreeMCMC(Options<Dynamics_>& o) : opts(o), dynamics(o)
    {
        std::cout<<"LikelihoodFreeMCMC(Options<Dynamics_>& o) Here."<<std::endl;  

        // Required options
        int K = opts.parallel_paths();
        int L = opts.path_length();
        int M = opts.extra_data_ratio();
        int N = opts.mcmc_trials();
        int c_dim = dynamics.parameter_dimension(o);
        
        std::cout<<"Here."<<std::endl;

        // The trajectories
        observed    = CoarsePathType( K, L, 2 );
        real        = PathType( K, L, M, 2 );
        x           = PathType( K, L, M, 2 );
        x_star      = PathType( K, L, M, 2 );

        // Parameters and chain
        c       = ParameterType(c_dim);
        c_star  = ParameterType(c_dim);
        real_c  = opts.parameters();
        chain   = ParameterChainType(N, c_dim);

        std::cout<<"Finished LikelihoodFreeMCMC(Options<Dynamics_>& o)"<<std::endl;
    }

    LikelihoodFreeMCMC()
    {
        std::cout<<"LikelihoodFreeMCMC() Here."<<std::endl;  
        // Default options.
        // opts = Options<Dynamics_>();
        // LikelihoodFreeMCMC(opts);
    }

    double log_acceptance_probability();
    double log_path_ratio();
    void propose_parameters(gsl_rng *r);
    void generate_true_trajectory(gsl_rng *r);
    void generate_observations(gsl_rng *r);
    void propose_trajectory(gsl_rng *r);
    void store_chain(int n);
    void accept();
    void finish();

private:
    Options<Dynamics_>      opts;   
    Dynamics_               dynamics;

    ParameterType           c;
    ParameterType           c_star;
    ParameterType           real_c;
    ParameterChainType      chain;
    
    CoarsePathType        observed;
    PathType              real;
    PathType              x;
    PathType              x_star;

};
      

template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_acceptance_probability()
{
    double log_total = 0;

    log_total += dynamics.log_transition_density_ratio(c_star, c);
    log_total += dynamics.log_prior_ratio(c_star, c);
    log_total += log_path_ratio();

    return log_total;
}
 
template<typename Dynamics_> 
void LikelihoodFreeMCMC<Dynamics_>::propose_parameters(gsl_rng *r)
    {
        int c_dim = dynamics.parameter_dimension(opts);
        for( size_t i=0; i<c_dim; ++i)
            c_star(i) = dynamics.sample_transition_density(r, c(i));
    }   

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    std::cout<<"begin LikelihoodFreeMCMC<Dynamics_>::generate_true_trajectory"<<std::endl;
    dynamics.generate_random_starts( r, real );
    std::cout<<"Done random starts..."<<std::endl;
    dynamics.trajectory(r, real_c, real );
    std::cout<<"leave LikelihoodFreeMCMC<Dynamics_>::generate_true_trajectory"<<std::endl;
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
    std::cout<<"begin LikelihoodFreeMCMC<Dynamics_>::generate_observations"<<std::endl;
    double variance = opts.observation_noise_variance();

    double random_noise_x;
    double random_noise_y;

    size_t K = opts.parallel_paths();
    size_t L = opts.path_length();

    for( size_t k=0; k<K; ++k )
    {
        for( size_t l=0; l<L; ++l )
        {   
            random_noise_y = gsl_ran_gaussian(r, variance );
            random_noise_x = gsl_ran_gaussian(r, variance );                

            observed(k, l, 0 ) = real(k, l, 0, 0 ) + random_noise_x;
            observed(k, l, 1 ) = real(k, l, 0, 1 ) + random_noise_y;
        }   
    }
    std::cout<<"leave LikelihoodFreeMCMC<Dynamics_>::generate_observations"<<std::endl;
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::propose_trajectory(gsl_rng *r) 
{
    dynamics.setup_observed_starts(r, observed, x_star );
    dynamics.trajectory(r, c_star, x_star ); 
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::store_chain(int n)
{
    int c_dim = dynamics.parameter_dimension(opts);
    for( size_t i=0; i<c_dim; ++i)
        chain(n, i) = c(i);

}
template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::accept() {
    x = x_star;
    c = c_star;
}
 
template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::finish() {
    dynamics.output_file_timeseries( chain );
}
          
template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_path_ratio()
{   
    double log_total = 0;
    log_total += dynamics.log_path_likelihood( x_star, c_star, observed );
    log_total -= dynamics.log_path_likelihood( x, c, observed );
    return log_total;
}
      
}



#endif
