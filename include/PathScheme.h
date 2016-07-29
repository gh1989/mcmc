#ifndef PATH_SCHEME_H
#define PATH_SCHEME_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"

namespace MCMC
{

class PathSchemeImp
{
    public:
        PathSchemeImp(Options &opts_) : opts(opts_) {}

    virtual void generate_true_trajectory( gsl_rng *r ) {}
    virtual void generate_observations( gsl_rng *r ) {}
    virtual void propose_trajectory( gsl_rng *r ) {}

    virtual double log_path_ratio() {return 1.0;}   
    virtual double log_transition_density_ratio() {return 1.0;}
    virtual double log_prior_ratio() {return 1.0;}

    virtual void propose_parameters( gsl_rng *r ) {}
    virtual void store_chain( int n ) {}
    virtual void accept() {}

    virtual void finish() {}

    protected:
        Options opts;

}; // end class PathSchemeImp


class PathScheme
{
    public:
    PathScheme(){};
    PathScheme(Options &opts_) {
        imp_ = new PathSchemeImp( opts_ );
    }
        
    virtual void generate_true_trajectory( gsl_rng *r ) {
        imp_->generate_true_trajectory(r);
    }

    virtual void generate_observations( gsl_rng *r ) {
        imp_->generate_observations(r);
    }

    virtual void propose_parameters( gsl_rng *r ) {
        imp_->propose_parameters(r);
    }

    virtual void propose_trajectory( gsl_rng *r ) {
        imp_->propose_trajectory(r);
    }

    virtual double log_path_ratio() {
        return imp_->log_path_ratio();
    }

    virtual double log_prior_ratio() {
        return imp_->log_prior_ratio();
    }

    virtual double log_transition_density_ratio() {
        return imp_->log_transition_density_ratio();
    }

    virtual void store_chain( int n ) {
        imp_->store_chain( n);
    }

    virtual void accept() {
        return imp_->accept();
    }

    virtual void finish() {
        imp_->finish();
    }        

    protected:
        PathSchemeImp *imp_; 

};  // end class PathScheme

}
#endif // PATH_SCHEME_H
