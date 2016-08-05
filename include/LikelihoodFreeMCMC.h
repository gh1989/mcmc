#ifndef LIKELIHOOD_MCMC_H
#define LIKELIHOOD_MCMC_H

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
class LikelihoodFreeMCMC
{
public:

    typedef typename Dynamics_::ParameterType ParameterType;
    typedef typename Dynamics_::ParameterChainType ParameterChainType;
    typedef typename Dynamics_::PathType PathType;
    typedef typename Dynamics_::CoarsePathType CoarsePathType;

    LikelihoodFreeMCMC(Options<Dynamics_>& o) : opts(o), dynamics(o)
    {
        // Required options
        K = opts.parallel_paths();
        L = opts.path_length();
        M = opts.extra_data_ratio();
        N = opts.mcmc_trials();
        c_dim = dynamics.parameter_dimension(o);

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
    }

    LikelihoodFreeMCMC()
    {
        // Default options.
        // opts = Options<Dynamics_>();
        // LikelihoodFreeMCMC(opts);
    }

    double log_acceptance_probability();

    double log_path_likelihood(     PathType &x, 
                                    ParameterType &c, 
                                    CoarsePathType &observed );

    

    void propose_parameters(gsl_rng *r);
    void propose_trajectory(gsl_rng *r);
    void propose(gsl_rng *r);
    void generate_true_trajectory(gsl_rng *r);
    void generate_observations(gsl_rng *r);

    void setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out );
    void generate_random_starts( gsl_rng *r, PathType &out );
    void trajectory( gsl_rng *r, ParameterType &c, PathType &out );
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
    
    size_t K;
    size_t L;
    size_t M;
    size_t N;
    size_t c_dim;

};
      

template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_acceptance_probability()
{
    double log_total = 0;

    log_total += dynamics.log_transition(c_star, c);
    log_total -= dynamics.log_transition(c, c_star);
    log_total += dynamics.log_prior(c_star);
    log_total -= dynamics.log_prior(c);
    log_total += log_path_likelihood( x_star, c_star, observed );
    log_total -= log_path_likelihood( x, c, observed );

    return log_total;
}
 
template<typename Dynamics_> 
void LikelihoodFreeMCMC<Dynamics_>::propose_parameters(gsl_rng *r)
{
    for( size_t i=0; i<c_dim; ++i)
        c_star(i) = dynamics.sample_transition_density(r, c(i));
}   

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::propose( gsl_rng *r )
{
    propose_parameters(r);
    propose_trajectory(r);
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    generate_random_starts( r, real );
    trajectory(r, real_c, real );
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
    double variance = opts.observation_noise_variance();

    double random_noise_x;
    double random_noise_y;

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
void LikelihoodFreeMCMC<Dynamics_>::setup_observed_starts( gsl_rng *r, 
                                                           CoarsePathType &y, 
                                                           PathType &out )
{
    double variance = opts.observation_noise_variance();
    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0); //+ gsl_ran_gaussian(r, variance);
        out(k, 0, 0, 1 ) = y(k, 0, 1); //+ gsl_ran_gaussian(r, variance);
    }
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::trajectory( gsl_rng *r, 
                                                ParameterType &c, 
                                                PathType &out )
{
    for( size_t k=0; k<K; ++k )
      for( size_t l=0; l<L-1; ++l )
      {
        for( size_t m=1; m<M; ++m )
        {
          dynamics.forward_sim(r, c, k, l, m, out);   
        }
        dynamics.forward_sim(r, c, k, l+1, 0, out);
      }
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::propose_trajectory(gsl_rng *r) 
{
    setup_observed_starts(r, observed, x_star );
    trajectory( r, c_star, x_star );
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::store_chain(int n)
{
    for( size_t i=0; i<c_dim; ++i)
        chain(n, i) = c(i);
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::accept() 
{
    x = x_star;   
    c = c_star;
}
 
template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::finish() 
{
    dynamics.output_file_timeseries( chain );
}

template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_path_likelihood(  PathType &x, 
                                                            ParameterType &c, 
                                                            CoarsePathType &observed )
{
    double log_total = 0;
    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L; ++l )
        {
            log_total += dynamics.log_path_likelihood( x, c, observed, k, l, 0 );
            for( size_t m=1; m<M; ++m )
            {
                log_total += dynamics.log_path_likelihood( x, c, observed, k, l, m );
            }
        }
    return log_total;
}
  
}

#endif
