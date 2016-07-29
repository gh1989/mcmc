#include "LikelihoodFreeScheme.h"

using namespace MCMC;

double LikelihoodFreeScheme::log_acceptance_probability()
{
    double log_total = 0;

    log_total += path_scheme.log_transition_density_ratio();
    log_total += path_scheme.log_prior_ratio();
    log_total += path_scheme.log_path_ratio();

    return log_total;
}
