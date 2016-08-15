#include "Dynamics.h"

using namespace MCMC;

double DynamicsBase::log_prior_sigma(double log_sigma)
{
    double proposal_sigma = _opts.parameter_proposal_diffusion_sigma();
    return -( 0.5 / (proposal_sigma*proposal_sigma) )*( pow(log_sigma, 2) );   
}

double DynamicsBase::log_transition_sigma(double log_sigma_star, double log_sigma)
{
    return 0;
}