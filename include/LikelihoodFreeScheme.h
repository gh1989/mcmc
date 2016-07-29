#ifndef LIKELIHOOD_FREE_SCHEME_H
#define LIKELIHOOD_FREE_SCHEME_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "PathScheme.h"

namespace MCMC
{

class LikelihoodFreeScheme
{
public:
    LikelihoodFreeScheme(const Options& opts_, const PathScheme& path_scheme_):opts(opts_), path_scheme(path_scheme_) {}

    LikelihoodFreeScheme(){}

    double log_acceptance_probability();

    void propose_parameters(gsl_rng *r){
        path_scheme.propose_parameters(r);
    }   

    void generate_true_trajectory(gsl_rng *r){ 
        path_scheme.generate_true_trajectory(r); 
    }

    void generate_observations(gsl_rng *r) { 
        path_scheme.generate_observations( r); 
    }
    
    void propose_trajectory(gsl_rng *r) { 
        path_scheme.propose_trajectory(r); 
    }

    void store_chain( int n ) {
        path_scheme.store_chain(n);
    }
   
    void accept() {
        path_scheme.accept();
    } 

    void finish() {
        path_scheme.finish();
    }
    
private:
    Options opts;   
    PathScheme path_scheme;
};
                           
}
#endif
