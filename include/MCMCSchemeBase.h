#ifndef MCMC_SCHEME_BASE_H
#define MCMC_SCHEME_BASE_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "PathScheme.h"

using namespace MCMC;

namespace MCMC
{

class MCMCSchemeBase
{

    public:
        MCMCSchemeBase();
        ~MCMCSchemeBase();

        virtual double log_acceptance_probability() = 0;
        virtual void propose_parameters(gsl_rng *r) = 0;
        virtual void generate_true_trajectory(gsl_rng *r) = 0;
        virtual void generate_observations(gsl_rng *r) = 0;
        virtual void propose_trajectory(gsl_rng *r) = 0;
        virtual void store_chain( int n ) = 0;
        virtual void accept() = 0;
        virtual void finish() = 0; 

};


}

#endif

