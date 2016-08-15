#ifndef SEQUENTIAL_LF_MCMC
#define SEQUENTIAL_LF_MCMC

/*
#include "Dynamics.h"
#include "MCMCBase.h"
#include "Options.h"


using namespace MCMC;

namespace MCMC
{

template<typename Dynamics_>
class SequentialLFMCMC
{
public:
    
    SequentialLFMCMC(Options& o) : (o)
    {
        // Samples from the SMC bridging distributions
        smc_trials = 100;
        trajectory_samples = SMCTrajectorySamplesType(smc_trials, K, L, M, 2);
        parameter_samples  = SMCParameterSamplesType(smc_trials, c_dim );
        sigma_samples      = SMCSigmaSamplesType(smc_trials);
    }

    SequentialLFMCMC() : MCMCBase<Dynamics_>() {}

    double log_acceptance_probability();
    void propose(gsl_rng *r);
    double log_path_likelihood( PathType &x, ParameterType &c, double sigma_, CoarsePathType &observed );
    void run_mcmc(gsl_rng *r, int l);


protected:
    size_t smc_trials;
    SMCTrajectorySamplesType trajectory_samples;
    SMCParameterSamplesType parameter_samples;
    SMCSigmaSamplesType sigma_samples;
    
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
            
            sigma_star = sigma_samples(random_int);
            
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
    
    double _sigma;
    double _sigma_star;
    
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
            _sigma = sigma_samples(i);
            
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
            _sigma = gsl_ran_lognormal( r, zeta_sigma, variance_sigma );
            for(size_t k=0; k<K; ++k)
            {
                for(size_t i=0; i<2; ++i)
                    _x(k, l, 0, i) = observed(k, l, i); // + gsl_ran_gaussian(r, variance);
            }
        }
        
        for(size_t j=0; j<c_dim; ++j)
            _c_star(j) = dynamics.sample_transition_density(r, _c(j));
        _sigma_star = gsl_ran_lognormal( r, log(_sigma), variance_sigma );
        
        for(size_t k=0; k<K; ++k)
        {
            _x_star(k, l, 0, 0) = _x(k, l, 0, 0);
            _x_star(k, l, 0, 1) = _x(k, l, 0, 1);
        
            for(size_t m=1; m<M; ++m)
            {
                dynamics.forward_sim(r, _c_star, _sigma_star, k, l, m, _x_star);
            }
            dynamics.forward_sim(r, _c_star, _sigma_star, k, l+1, 0, _x_star);
        }
        log_u = log( gsl_rng_uniform(r) );
        
        log_a = 0;
        for(size_t k=0; k<K; ++k)
        {
            log_a += dynamics.log_path_likelihood( _x_star, _c_star, _sigma_star, observed, k, l+1, 0 );
            log_a -= dynamics.log_path_likelihood( _x, _c, _sigma, observed, k, l+1, 0 );
           
        }
        if( log_u < log_a )
        {
            _x = _x_star;
            _c = _c_star;
            _sigma = _sigma_star;
        }
        
        for( size_t k=0; k<K; ++k)
        for( size_t m=0; m<M; ++m)
        for( size_t j=0; j<2; ++j)
            trajectory_samples(i, k, l, m, j) = _x(k, l, m, j);
        
        for( size_t j=0; j<c_dim; ++j)
            parameter_samples(i, j) = _c(j);
        sigma_samples(i) = _sigma;
        
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
    log_total += log_path_likelihood( x_star, c_star, sigma_star, observed );
    log_total -= log_path_likelihood( x, c, sigma, observed );

    return log_total;
}

template<typename Dynamics_>
double SequentialLFMCMC<Dynamics_>::log_path_likelihood(  PathType &x, 
                                                          ParameterType &c, 
                                                          double sigma_,
                                                          CoarsePathType &observed )
{
    double log_total = 0;
    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L; ++l )
        {
            log_total += dynamics.log_path_likelihood( x, c, sigma_, observed, k, l, 0 );
            for( size_t m=1; m<M; ++m )
            {
                log_total += dynamics.log_path_likelihood( x, c, sigma_, observed, k, l, m );
            }
        }
    return log_total;
}

}
*/
#endif
