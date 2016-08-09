#ifndef MCMC_BASE_H
#define MCMC_BASE_H

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

using namespace MCMC;

namespace MCMC
{
 
template< class Dynamics_ > 
class MCMCBase
{

public:
    typedef typename Dynamics_::ParameterType ParameterType;
    typedef typename Dynamics_::ParameterChainType ParameterChainType;
    typedef typename Dynamics_::PathType PathType;
    typedef typename Dynamics_::CoarsePathType CoarsePathType;
    typedef Tensor<double, 1> SigmaChainType;

    MCMCBase(Options &o) :  opts(o), dynamics(o)
    {        
        // MCMC options
        K = o.parallel_paths();
        L = o.path_length();
        M = o.extra_data_ratio();
        N = o.mcmc_trials();

        // The trajectories
        observed    = CoarsePathType( K, L, 2 );
        real        = PathType( K, L, M, 2 );
        x           = PathType( K, L, M, 2 );
        x_star      = PathType( K, L, M, 2 );     
        
     
        
        // Parameters
        c_dim               = dynamics.parameter_dimension(o);
        c                   = ParameterType(c_dim);
        
           for(size_t i=0; i<c_dim; ++i)
            std::cout<< "c " << i << " "<< c(i) << std::endl;
        
        c_star              = ParameterType(c_dim);
        real_c              = dynamics.default_parameters(o);
        parameter_chain     = ParameterChainType(N, c_dim);
        
        // Parameter: diffusion coefficient
        real_sigma          = o.real_sigma();
        sigma               = o.start_sigma();
        zeta_sigma          = o.zeta_sigma();
        variance_sigma      = o.variance_sigma();  
        sigma_chain         = SigmaChainType(N);
    }
    
    MCMCBase(){}
    ~MCMCBase(){}
        
    void generate_true_trajectory(gsl_rng *r);
    void generate_random_starts( gsl_rng *r, PathType &out );
    void generate_observations(gsl_rng *r);
    void trajectory( gsl_rng *r, ParameterType &c, double sigma_, PathType &out );
    void setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out );
    void store_chain(int n);
    void accept();
    void finish();
    
    virtual double log_acceptance_probability() = 0;
    virtual void propose(gsl_rng *r) = 0;
    
protected:
    
    Options                 opts;   
    Dynamics_               dynamics;
    
    size_t K;
    size_t L;
    size_t M;
    size_t N;

    CoarsePathType          observed;
    PathType                real;
    PathType                x;
    PathType                x_star;
    
    size_t                  c_dim;
    ParameterType           c;
    ParameterType           c_star;
    ParameterType           real_c;
    ParameterChainType      parameter_chain;
    
    double                  real_sigma;
    double                  zeta_sigma;
    double                  sigma_star;
    double                  sigma; 
    double                  variance_sigma;    
    SigmaChainType          sigma_chain; 
}; 
    
}


template<typename Dynamics_>
void MCMCBase<Dynamics_>::store_chain(int n)
{
    for( size_t i=0; i<c_dim; ++i)
        parameter_chain(n, i) = c(i);
    sigma_chain(n) = sigma;
}

template<typename Dynamics_>
void MCMCBase<Dynamics_>::accept() 
{
    x = x_star;   
    c = c_star;
    sigma = sigma_star;
}
 
template<typename Dynamics_>
void MCMCBase<Dynamics_>::finish() 
{
    dynamics.output_file_timeseries( parameter_chain );
}


template<typename Dynamics_>
void MCMCBase<Dynamics_>::setup_observed_starts(  gsl_rng *r, 
                                                  CoarsePathType &y, 
                                                  PathType &out )
{
    //double variance = opts.observation_noise_variance();
    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0); //+ gsl_ran_gaussian(r, variance);
        out(k, 0, 0, 1 ) = y(k, 0, 1); //+ gsl_ran_gaussian(r, variance);
    }
}

template<typename Dynamics_>
void MCMCBase<Dynamics_>::trajectory( gsl_rng *r, 
                                      ParameterType &c,
                                      double sigma_,
                                      PathType &out )
{
    for( size_t k=0; k<K; ++k )
      for( size_t l=0; l<L-1; ++l )
      {
        for( size_t m=1; m<M; ++m )
        {
          dynamics.forward_sim(r, c, sigma_, k, l, m, out);   
        }
        dynamics.forward_sim(r, c, sigma_, k, l+1, 0, out);
        std::cout<<"out "<<l<< out(0,l+1,0,0)<<" , "<< out(0,l+1,0,1)<<std::endl;
      }
}

template<typename Dynamics_>
void MCMCBase<Dynamics_>::generate_observations(gsl_rng *r)
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
            
            std::cout<<"observed" << l << " "<< observed(k, l, 0 ) << ","<< observed(k, l, 1 ) << std::endl;
        }   
    }
}

template<typename Dynamics_>
void MCMCBase<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

template<typename Dynamics_>
void MCMCBase<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    generate_random_starts( r, real );
    trajectory(r, real_c, real_sigma, real );
}

#endif