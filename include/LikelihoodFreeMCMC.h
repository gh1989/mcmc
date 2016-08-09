#ifndef LIKELIHOOD_MCMC_H
#define LIKELIHOOD_MCMC_H

#include "Dynamics.h"
#include "MCMCBase.h"
#include "Options.h"

using namespace MCMC;

namespace MCMC
{

template<typename Dynamics_>
class LikelihoodFreeMCMC : public MCMCBase<Dynamics_>
{
public:
    typedef typename Dynamics_::ParameterType ParameterType;
    typedef typename Dynamics_::ParameterChainType ParameterChainType;
    typedef typename Dynamics_::PathType PathType;
    typedef typename Dynamics_::CoarsePathType CoarsePathType;
    typedef Tensor<double, 1> SigmaChainType; 
    
    using MCMCBase<Dynamics_>::opts;
    using MCMCBase<Dynamics_>::dynamics;
       
    using MCMCBase<Dynamics_>::K;
    using MCMCBase<Dynamics_>::M;
    using MCMCBase<Dynamics_>::L;
    using MCMCBase<Dynamics_>::N;
    
    using MCMCBase<Dynamics_>::observed;
    using MCMCBase<Dynamics_>::real;
    using MCMCBase<Dynamics_>::x;
    using MCMCBase<Dynamics_>::x_star;
    
    using MCMCBase<Dynamics_>::c_dim;
    using MCMCBase<Dynamics_>::c;
    using MCMCBase<Dynamics_>::c_star;
    using MCMCBase<Dynamics_>::real_c;
    using MCMCBase<Dynamics_>::parameter_chain;
    
    using MCMCBase<Dynamics_>::real_sigma;
    using MCMCBase<Dynamics_>::zeta_sigma;
    using MCMCBase<Dynamics_>::sigma_star;
    using MCMCBase<Dynamics_>::sigma;
    using MCMCBase<Dynamics_>::variance_sigma;
    using MCMCBase<Dynamics_>::sigma_chain;
   
    using MCMCBase<Dynamics_>::setup_observed_starts;
    using MCMCBase<Dynamics_>::trajectory;
   
    LikelihoodFreeMCMC(Options& o) : MCMCBase<Dynamics_>(o) {}
    LikelihoodFreeMCMC() : MCMCBase<Dynamics_>() {}
    
    /* virtual *////////////////////////
    void propose(gsl_rng *r); 
    double log_acceptance_probability();
    ////////////////////////////////////
    
    double log_path_likelihood(PathType &path, ParameterType &params, double sigma_ );
};
 
} 

template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_acceptance_probability()
{
    double log_total = 0;
    
    log_total += dynamics.log_transition(c_star, c);
    log_total -= dynamics.log_transition(c, c_star);
    log_total += log( gsl_ran_lognormal_pdf(sigma_star, zeta_sigma, variance_sigma) );
    log_total -= log( gsl_ran_lognormal_pdf(sigma, zeta_sigma, variance_sigma) );
    log_total += dynamics.log_prior(c_star);
    log_total -= dynamics.log_prior(c);
    log_total += log_path_likelihood( x_star, c_star, sigma_star );
    log_total -= log_path_likelihood( x, c, sigma );

    return log_total;
}

template<typename Dynamics_>
void LikelihoodFreeMCMC<Dynamics_>::propose( gsl_rng *r )
{
    for( size_t i=0; i<c_dim; ++i)
        c_star(i) = dynamics.sample_transition_density(r, c(i));
    sigma_star = gsl_ran_lognormal( r, log(sigma), variance_sigma );
    for(size_t l=0; l<L; ++l)
        std::cout<< "x_star " <<l << " "<< x_star(0,l,0,0) <<", "<< x_star(0,l,0,1) << std::endl;
    setup_observed_starts(r, observed, x_star );
    for(size_t l=0; l<L; ++l)
        std::cout<< "x_star " <<l << " "<< x_star(0,l,0,0) <<", "<< x_star(0,l,0,1) << std::endl;
    trajectory( r, c_star, sigma_star, x_star );
    for(size_t l=0; l<L; ++l)
        std::cout<< "x_star " <<l << " "<< x_star(0,l,0,0) <<", "<< x_star(0,l,0,1) << std::endl;
}

template<typename Dynamics_>
double LikelihoodFreeMCMC<Dynamics_>::log_path_likelihood(  PathType &path, 
                                                            ParameterType &c_, 
                                                            double sigma_ )
{
    double log_total = 0;
    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L; ++l )
        {
            log_total += dynamics.log_path_likelihood( path, c_, sigma_, observed, k, l, 0 );
            for( size_t m=1; m<M; ++m )
            {
                log_total += dynamics.log_path_likelihood( path, c_, sigma_, observed, k, l, m );
            }
        }
    std::cout<<"LikelihoodFreeMCMC@log_path_likelihood. log_total: "<< log_total << std::endl;
    return log_total;
}

#endif
