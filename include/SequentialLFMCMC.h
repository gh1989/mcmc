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
    typedef typename Dynamics_::ParameterPointType ParameterPointType;
    
    typedef Tensor<double, 5> SMCTrajectorySamplesType;
    typedef Tensor<ParameterPointType, 2> SMCParameterSamplesType;
    
    SequentialLFMCMC(Options<Dynamics_>& o) : opts(o), dynamics(o)
    {
        // Required options
        K = opts.parallel_paths();
        L = opts.path_length();
        M = opts.extra_data_ratio();
        N = opts.mcmc_trials();
        smc_trials = N;
        c_dim = dynamics.parameter_dimension(o);

        // The trajectories
        observed    = CoarsePathType( K, L, 2 );
        real        = PathType( K, L, M, 2 );
        x           = PathType( K, L, M, 2 );
        x_star      = PathType( K, L, M, 2 );

        // Samples from the SMC bridging distributions
        trajectory_samples = SMCTrajectorySamplesType(smc_trials, K, L, M, 2);
        parameter_samples = SMCParameterSamplesType(smc_trials, c_dim );
        
        // Parameters and chain
        c       = ParameterType(c_dim);
        c_star  = ParameterType(c_dim);
        real_c  = opts.parameters();
        chain   = ParameterChainType(N, c_dim);
    }

    SequentialLFMCMC()
    {
        // Default options.
        // opts = Options<Dynamics_>();
        // SequentialLFMCMC(opts);
    }

    double log_acceptance_probability();

    double log_path_likelihood(     PathType &x, 
                                    ParameterType &c, 
                                    CoarsePathType &observed );

    

    void propose(gsl_rng *r);
    void run_mcmc(gsl_rng *r, int l);
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
    
    SMCTrajectorySamplesType    trajectory_samples;
    SMCParameterSamplesType     parameter_samples;
    
    size_t K;
    size_t L;
    size_t M;
    size_t N;
    size_t smc_trials;
    size_t c_dim;

};
      
template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::propose( gsl_rng *r )
{
    int random_int;
    random_int = gsl_rng_uniform_int(r, smc_trials);
    
    for(size_t l=0; l<L-1; ++l)
        run_mcmc(r, l);
            
    for( size_t k=0; k<K; ++k)
    {      
        for(size_t l=0; l<L; ++l)
        {
            // For a given k the c_star will be?
            for(size_t i=0; i<c_dim; ++i)
                c_star(i) = parameter_samples(random_int, i);
            
            for(size_t m=0; m<M; ++m)
                for(size_t i=0; i<2; ++i)
                    x_star(k,l,m,i) = trajectory_samples(random_int, k, l, m, i);
                
        }
    }
   
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::run_mcmc(gsl_rng *r, int l)
{   
    int random_int;
    
    ParameterType _c(c_dim);
    ParameterType _c_star(c_dim);
    
    // K*L*M*2 sized Tensor but only using one time slot. Compatibility.
    PathType _x(K, L, M, 2);
    PathType _x_star(K, L, M, 2);

    double log_u;
    double log_a;
    
    for( size_t i=0; i<smc_trials; ++i)
    {
        if( l>0 )
        {
            // Initialise from (t-1)
            random_int = gsl_rng_uniform_int(r, smc_trials);
            for(size_t i=0; i<c_dim; ++i)
                _c(i) = parameter_samples(random_int, i);
            for(size_t k=0; k<K; ++k)
            {
                for(size_t i=0; i<2; ++i)
                    _x(k, l, 0, i) = trajectory_samples(random_int, k, l, 0, i);
            }
        }
        else
        {
            // Initialise from the prior.
            for(size_t i=0; i<c_dim; ++i)
                _c(i) = dynamics.sample_transition_density(r, 0);
            for(size_t k=0; k<K; ++k)
            {
                for(size_t i=0; i<2; ++i)
                    _x(k, l, 0, i) = observed(k, l, i); // + gsl_ran_gaussian(r, variance);
            }
        }
        
        for(size_t j=0; j<c_dim; ++j)
            _c_star(j) = dynamics.sample_transition_density(r, _c(j));
        
        for(size_t k=0; k<K; ++k)
        {
            _x_star(k, l, 0, 0) = _x(k, l, 0, 0);
            _x_star(k, l, 0, 1) = _x(k, l, 0, 1);
        
            for(size_t m=1; m<M; ++m)
            {
                dynamics.forward_sim(r, _c_star, k, l, m, _x_star);
            }
            dynamics.forward_sim(r, _c_star, k, l+1, 0, _x_star);
        }
        log_u = log( gsl_rng_uniform(r) );
        
        log_a = 0;
        for(size_t k=0; k<K; ++k)
        {
            log_a += dynamics.log_path_likelihood( _x_star, _c_star, observed, k, l+1, 0 );
            log_a -= dynamics.log_path_likelihood( _x, _c, observed, k, l+1, 0 );
           
        }
        if( log_u < log_a )
        {

            _x = _x_star;
            _c = _c_star;
        }
        
        for( size_t k=0; k<K; ++k)
        for( size_t m=0; m<M; ++m)
        for( size_t j=0; j<2; ++j)
            trajectory_samples(i, k, l, m, j) = _x(k, l, m, j);
        
        for( size_t j=0; j<c_dim; ++j)
            parameter_samples(i, j) = _c(j);
        
    }
    
}
      
template<typename Dynamics_>
double SequentialLFMCMC<Dynamics_>::log_acceptance_probability()
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
void SequentialLFMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    generate_random_starts( r, real );
    trajectory(r, real_c, real );
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::generate_observations(gsl_rng *r)
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
void SequentialLFMCMC<Dynamics_>::setup_observed_starts( gsl_rng *r, 
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
void SequentialLFMCMC<Dynamics_>::trajectory( gsl_rng *r, 
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
void SequentialLFMCMC<Dynamics_>::store_chain(int n)
{
    for( size_t i=0; i<c_dim; ++i)
        chain(n, i) = c(i);
}

template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::accept() 
{
    x = x_star;   
    c = c_star;
}
 
template<typename Dynamics_>
void SequentialLFMCMC<Dynamics_>::finish() 
{
    dynamics.output_file_timeseries( chain );
}

template<typename Dynamics_>
double SequentialLFMCMC<Dynamics_>::log_path_likelihood(  PathType &x, 
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
/*
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
