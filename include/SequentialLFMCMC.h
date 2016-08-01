#ifndef SEQUENTIAL_LF_MCMC
#define SEQUENTIAL_LF_MCMC

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "Dynamics.h"

using namespace MCMC;

namespace MCMC
{

template< typename Dynamics_ >
class SequentialLFMCMC
{
};

}
/*
void SequentialLFMCMC::sequential_step()
{
    if( t == 0 )
    {
        // Sample from the prior
        parameters = dynamics.params().sample_parameter_prior();
        xt = dynamics.sample_path_prior();
    }
    else
    {
        // Sample from the (large) sample from p( theta, x_t | Y_t )
        parameters = sample_parameters();
        xt = sample_path();
    }
    
    parameters = dynamics.params().propose( parameters );
    path = dynamics.forward_simulate( xt );  
    xtplus1 = path(M-1);

    for( size_t i=0; i<N; ++i )
    {
        u = gsl_ran_uniform(r);
        // Is this caled the marginal? p(y(t+1) | x*(t+1), c* )
        A = 0;
        A += dynamics.log_marginal( observed(t+1), xtplus1_star, parameters_star );
        A -= dynamics.log_marginal( observed(t+1), xtplus1,      parameters );
    }

}
*/
#endif
