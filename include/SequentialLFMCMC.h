#ifndef SEQUENTIAL_LF_MCMC
#define SEQUENTIAL_LF_MCMC

#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

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
class SequentialLFMCMC
{
public:

    typedef typename Dynamics_::ParameterType ParameterType;
    typedef typename Dynamics_::ParameterChainType ParameterChainType;
    typedef typename Dynamics_::PathType PathType;
    typedef typename Dynamics_::CoarsePathType CoarsePathType;

    SequentialLFMCMC(Options<Dynamics_>& o) : opts(o), dynamics(o)
    {
        // Required options
        int K = opts.parallel_paths();
        int L = opts.path_length();
        int M = opts.extra_data_ratio();
        int N = opts.mcmc_trials();
        int c_dim = dynamics.parameter_dimension(o);

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

        // SMC sub_parameters and samples
        // smc_samples_c = Tensor< ComplexType, >( , L, M ,c_dim);
        // smc_samples_x = Tensor< double, >(K, L, M ,2)
    }

    SequentialLFMCMC()
    {
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
    size_t _smc_mcmc_trials = 100;

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
double SequentialLFMCMC<Dynamics_>::log_acceptance_probability()
{
    double log_total = 0;

    log_total += dynamics.log_transition_density_ratio(c_star, c);
    log_total += dynamics.log_prior_ratio(c_star, c);
    log_total += log_path_ratio();

    return log_total;
}
 
template<typename Dynamics_> 
void SequentialLFMCMC<Dynamics_>::propose_parameters(gsl_rng *r)
    {
        int c_dim = dynamics.parameter_dimension(opts);
        for( size_t i=0; i<c_dim; ++i)
            c_star(i) = dynamics.sample_transition_density(r, c(i));
    }   

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    dynamics.generate_random_starts( r, real );
    dynamics.trajectory(r, real_c, real );

}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
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
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::propose_trajectory(gsl_rng *r) 
{
    dynamics.setup_observed_starts(r, observed, x_star );
    dynamics.trajectory(r, c_star, x_star ); 
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::store_chain(int n)
{
    int c_dim = dynamics.parameter_dimension(opts);
    for( size_t i=0; i<c_dim; ++i)
        chain(n, i) = c(i);

}
template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::accept() 
{
    int K = opts.parallel_paths();
    int L = opts.path_length();
    int M = opts.extra_data_ratio();
    
    for( size_t k=0; k<K; ++k)
        for( size_t l=0; l<L; ++l)
            for( size_t m=0; m<M; ++m)
                for( size_t i=0; i<2; ++i)
                    x(k,l,m,i) = x_star(k,l,m,i);
    
    c = c_star;
}
 
template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::finish() 
{
    dynamics.output_file_timeseries( chain );
}
          
template<typename Dynamics_>
double SequentialLFMCMC<Dynamics_>::log_path_ratio()
{   
    double log_total = 0;
    log_total += dynamics.log_path_likelihood( x_star, c_star, observed );
    log_total -= dynamics.log_path_likelihood( x, c, observed );
      
    return log_total;
}
 /*
template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::bridge_distributions( gsl_rng *r, size_t t )
{
   
    int K = opts.parallel_paths();
    int L = opts.path_length();
    int M = opts.extra_data_ratio();
    int N = opts.mcmc_trials();
    int c_dim = dynamics.parameter_dimension();

    double log_acceptance_ratio;
    double observation_variance = opts.observation_variance();
    double parameter_variance = opts.paramteter_variance(); 

    ParameterType cc( c_dim );
    PathType xt( K, L, M, 2 );

    xtplus1;
    xtplus1_star;
    cc;
    cc_star;

    if( t==0 ) {
        cc = gsl_rng_gaussian(r, parameter_variance);
        xt = gsl_rng_guassian(r, observation_variance) + observed(0);
    }
    else {
        cc = smc_samples_c( t, random_index );
        xt = smc_samples_x( t, random_index );
    }

    // Run an MCMC targetting 
    for( size_t i=0; i<_smc_mcmc_trials; ++i )
    {
        // Perturb
        for( size_t i=0; i<c_dim; ++i)
            cc_star(i) = dynamics.sample_transition_density(r, c(i));

        // Forward simulate
        dynamics.trajectory( r, xt, inbetween_path );
        xx_star = inbetween_path(M-1);
        
        log_acceptance_ratio = 0
        log_total += dynamics.log_path_likelihood( xtplus1_star, cc_star, observed );
        log_total -= dynamics.log_path_likelihood( xtplus1, cc, observed );       

        log_u = log( gsl_ran_uniform(r) );
        if( log_u < log_acceptance_ratio ) {
            xx = xx_star;
            cc = cc_star;
        }       

        for(size_t j=0; j<c_dim; ++j)
            smc_samples_c( t+1, i, j) = cc(j);
        for(size_t j=0; j<2; ++j )
            smc_samples_x( t+1, i, j) = xx(j);
        
    }
*/
   
      
}

#endif
