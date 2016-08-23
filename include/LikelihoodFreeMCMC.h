#ifndef LIKELIHOOD_MCMC_H
#define LIKELIHOOD_MCMC_H

#ifndef M_PI
    #define M_PI 3.14159265359
#endif
#include <chrono>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Options.h"

using namespace Eigen;
using namespace MCMC;

namespace MCMC
{

template<class Dynamics_>
class LikelihoodFreeMCMC
{
public:

    using ParameterType         = typename Dynamics_::ParameterType;
    using ParameterChainType    = typename Dynamics_::ParameterChainType;
    using PathType              = typename Dynamics_::PathType;
    using CoarsePathType        = typename Dynamics_::CoarsePathType;
    using SigmaChainType        = typename Dynamics_::SigmaChainType; 
   
    LikelihoodFreeMCMC(Options& o) : _dynamics(o)
    {
        std::cout<< this << "LikelihoodFreeMCMC(Options& o)" << std::endl;
        _opts.print_options(std::cout);
        setup_from_options(o);
    }
    
    // Do not instantiate without options.
    LikelihoodFreeMCMC()=delete;
    
    LikelihoodFreeMCMC& operator=(LikelihoodFreeMCMC&& other)
    {
        std::cout<<"LikelihoodFreeMCMC move assignment operator= called in " << this <<std::endl;
        _dynamics = other.dynamics();
        setup_from_options( other.opts() );
        return *this;
    }    
    
    LikelihoodFreeMCMC(LikelihoodFreeMCMC&& other) : _dynamics(other.dynamics())
    {
        std::cout<<"LikelihoodFreeMCMC move constructor called in " << this <<std::endl;
        setup_from_options( other.opts() );
    }
    
    LikelihoodFreeMCMC& operator=(LikelihoodFreeMCMC& other)
    {
        std::cout<<"LikelihoodFreeMCMC move assignment operator= called in " << this <<std::endl;
        _dynamics = other.dynamics();
        setup_from_options( other.opts() );
        return *this;
    }    
    
    LikelihoodFreeMCMC(LikelihoodFreeMCMC& other) : _dynamics(other.dynamics())
    {
        std::cout<<"LikelihoodFreeMCMC copy constructor called in " << this <<std::endl;
        setup_from_options( other.opts() );
    }
    
    ~LikelihoodFreeMCMC()
    {
        std::cout<<"~LikelihoodFreeMCMC() called in " << this <<std::endl;
    };
    
    void accept();
    void finish();
    void generate_observations( gsl_rng *r);
    void generate_true_trajectory(gsl_rng *r);
    void propose(gsl_rng *r); 
    void store_chain(int n);
    void setup_from_options(Options &o);
    
    double log_acceptance_probability();
    double log_path_likelihood(PathType &path, ParameterType &params, double sigma_ );
    
    Options& opts(){ return _opts; }
    Dynamics_& dynamics(){return _dynamics;}
    
    ParameterType& current_drift(){ return c; }
    double current_log_sigma() { return log_sigma; }
    
private:
    void setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out );
    void generate_random_starts( gsl_rng *r, PathType &out );
    void trajectory( gsl_rng *r, ParameterType &c, double sigma_, PathType &out );
    
    Options                 _opts;   
    Dynamics_               _dynamics;
    
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
    ParameterType           c_average;
    
    double                  log_real_sigma;
    double                  log_sigma_star;
    double                  log_sigma;  
    SigmaChainType          log_sigma_chain; 
    double                  log_sigma_average;
    
    bool    infer_drift_parameters;
    bool    infer_diffusion_parameters;
    
    std::ofstream mcmc_file;
    std::string filename;
};
 
} 

template<class Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_acceptance_probability()
{
    double log_total = 0;
       
    // b
    if( infer_drift_parameters )
    {      
        log_total += _dynamics.log_transition(c, c_star);       
        log_total -= _dynamics.log_transition(c_star, c);
        log_total += _dynamics.log_prior(c_star);
        log_total -= _dynamics.log_prior(c);

    }

    // Sigma
    if( infer_diffusion_parameters )
    {
        log_total += _dynamics.log_transition_sigma(log_sigma, log_sigma_star);
        log_total -= _dynamics.log_transition_sigma(log_sigma_star, log_sigma);
        log_total += _dynamics.log_prior_sigma(log_sigma_star);
        log_total -= _dynamics.log_prior_sigma(log_sigma);
    }
    
    log_total += log_path_likelihood( x_star, c_star, log_sigma_star );
    log_total -= log_path_likelihood( x, c, log_sigma );

    return log_total;
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::propose( gsl_rng *r )
{   
    if( infer_drift_parameters )
    {
        //for( size_t i=0; i<c_dim; ++i)
        c_star(c_dim-1) = _dynamics.sample_transition_density(r, c(c_dim-1));
    }
    else
        c_star = c;
    
    double log_sigma_proposal_standard_deviation = _opts.parameter_proposal_diffusion_sigma();
    
    if( infer_diffusion_parameters )
        log_sigma_star = gsl_ran_gaussian(r, log_sigma_proposal_standard_deviation ) + log_sigma;
    else
        log_sigma_star = log_sigma;
    
    setup_observed_starts(r, observed, x_star );
    trajectory( r, c_star, log_sigma_star, x_star );
}

template<class Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_path_likelihood(  PathType &path, ParameterType &c_,  double log_sigma_ )
{
    double log_total = 0;
    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L; ++l )
        {
            log_total += _dynamics.log_path_likelihood( path, c_, log_sigma_, observed, k, l, 0 );
            if( l < L-1 )
            {
                for( size_t m=1; m<M; ++m )
                    log_total += _dynamics.log_path_likelihood( path, c_, log_sigma_, observed, k, l, m );
            }
        }
    return log_total;
}


template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::store_chain(int n)
{
    if (_opts.store_time_series())
    {
        for( size_t i=0; i<c_dim; ++i)
            parameter_chain(n, i) = c(i);
        log_sigma_chain(n) = log_sigma;
    }
    else
    {
        //c_average += c;
        log_sigma_average += log_sigma;
    }
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::accept() 
{
    x = x_star;   
    c = c_star;
    log_sigma = log_sigma_star;
}
 
template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::finish() 
{
    // experiment options
    mcmc_file.open( filename );
    
    if (_opts.store_time_series())
    {

    _opts.print_header( mcmc_file );
    _dynamics.output_file_timeseries( parameter_chain, log_sigma_chain, mcmc_file );
    
    }
    else
    {
        log_sigma_average /= N;
    
        //for( size_t i=0; i<c_dim; ++i)
        //    c_average(i) = c_average(i)/double(N);
    
        mcmc_file << "log_sigma_average" << std::endl;
        mcmc_file << log_sigma_average << std::endl;
    
        //mcmc_file << "c_average" << std::endl;
        //mcmc_file << c_average << std::endl;
    }
    
    mcmc_file.close();
}


template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::setup_observed_starts(  gsl_rng *r, 
                                                            CoarsePathType &y, 
                                                            PathType &out )
{
    double sigma = _opts.observation_noise_sigma();
    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0);
        out(k, 0, 0, 1 ) = y(k, 0, 1);
    }
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::trajectory(   gsl_rng *r, 
                                                  ParameterType &c_,
                                                  double log_sigma_,
                                                  PathType &out )
{
    for( size_t k=0; k<K; ++k )
      for( size_t l=0; l<L-1; ++l )
      {
        for( size_t m=1; m<M; ++m )
        {
          _dynamics.forward_sim(r, c_, log_sigma_, observed, k, l, m, out);   
        }
        _dynamics.forward_sim(r, c_, log_sigma_, observed, k, l+1, 0, out);
      }
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
    double observation_sigma = _opts.observation_noise_sigma();
    //std::cout << "LikelihoodFreeMCMC thinks that observation_noise_sigma is: " << observation_sigma << std::endl;
    double random_noise_x;
    double random_noise_y;

    for( size_t k=0; k<K; ++k )
    {
        for( size_t l=0; l<L; ++l )
        {   
            random_noise_y = gsl_ran_gaussian(r, observation_sigma );
            random_noise_x = gsl_ran_gaussian(r, observation_sigma );                

            observed(k, l, 0 ) = real(k, l, 0, 0 ) + random_noise_x;           
            observed(k, l, 1 ) = real(k, l, 0, 1 ) + random_noise_y;
        }   
    }
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    generate_random_starts( r, real );
    trajectory(r, real_c, log_real_sigma, real );
}

template<class Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::setup_from_options(Options &o)
{
    _opts = o;

    // The number of parallel trajectories per iteration.
    K = o.parallel_paths();
    
    // The length of each trajectory path in timesteps.
    L = o.path_length();
    
    // The extra data in between observed data ratio.
    M = o.extra_data_ratio();
    
    // The number of iterations for the MCMC experiment.
    N = o.mcmc_trials();

    // The observations K*L*2 Tensor of doubles
    observed    = CoarsePathType( K, L, 2 );
    
    // The actual path
    real        = PathType( K, L, M, 2 );
    
    // The proposed and current path for the MCMC simulation
    x           = PathType( K, L, M, 2 );
    x_star      = PathType( K, L, M, 2 );    
    
    // We may want to infer these separately
    infer_drift_parameters = o.infer_drift_parameters();
    infer_diffusion_parameters = o.infer_diffusion_parameters();
    
    // The trajectories: dX_t = b(X_t) dt + sqrt(2*D(X_t)) dW_t
    // b(X_t) : related parameters
    
    // The dimension of the b(X_t) parameters
    c_dim               = _dynamics.parameter_dimension(o);
    
    // The current parameters for b(X_t)
    c                   = ParameterType(c_dim);
    // The proposed parameters for b(X_t)
    c_star              = ParameterType(c_dim);
    
    // The genuine parameters that need to be inferred
    real_c              = _dynamics.default_parameters(o);
    
    if( infer_drift_parameters )
        for(size_t i=0; i<c_dim; ++i)  c(i) = 0;
    else
        c = real_c;
            
    // The chain of b(X_t) parameters
    parameter_chain     = ParameterChainType(N, c_dim);
        
    // D(X_t) : in this case just the diffusion coefficient, distributed log normally so we
    // propose log_sigma with Gaussian prior
    log_real_sigma              = o.log_real_sigma();
    
    if( infer_diffusion_parameters )
    {
        std::cout<<"Inferring diffusion parameters."<<std::endl;
        log_sigma               = o.log_start_sigma();
    }
    else
        log_sigma = log_real_sigma;
    
    log_sigma_chain         = SigmaChainType(N);    
    
    unsigned long int milliseconds_since_epoch = 
    std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
        
    filename = "output/";
    filename += _opts.output_subfolder();
    filename += "/LFTimeSeries_";
    filename += Dynamics_::dynamics_string();
    filename += std::to_string( milliseconds_since_epoch );
    filename += ".txt";

}

#endif
