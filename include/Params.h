#ifndef PARAMS_H
#define PARAMS_H

#include "FourierSeries.h"
#include "Options.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace MCMC
{

class ParamImp {
    public:
        ParamImp(Options &opts_) : opts(opts_) {}
        ParamImp(){}

    virtual void propose_parameters(gsl_rng *r) {}
    virtual double log_transition_density_ratio() { return 1.0; }
    virtual double log_prior_ratio() { return 1.0; }
    virtual void accept() {}
    virtual void store_chain(int n) {}
    virtual void output_file_timeseries() {}

    // Work around...
    virtual FourierSeries potential(){ return FourierSeries(1); }
    virtual FourierSeries potential_star(){ return FourierSeries(1); }
    virtual double c(){ return 1.0; }
    virtual double c_star(){ return 1.0; }

    protected:
        Options opts;
};

class Params {
  public:
    Params(){}
    Params(Options &opts) {
        imp_ = new ParamImp(opts);
    }
    virtual void propose_parameters(gsl_rng *r) {
        imp_->propose_parameters(r);
    }    

    virtual void output_file_timeseries(){
        imp_->output_file_timeseries();
    }

    virtual double log_transition_density_ratio() {
        return imp_->log_transition_density_ratio();
    }

    virtual double log_prior_ratio() {
        return imp_->log_prior_ratio();
    }

    virtual void accept() {
        imp_->accept();
    }

    virtual void store_chain(int n) {
        imp_->store_chain(n);
    }

    // Work around...
    virtual FourierSeries potential(){
         return imp_->potential(); 
    }
    
    virtual FourierSeries potential_star(){ 
        return imp_->potential_star(); 
    }

    virtual double c(){ return imp_->c(); }
    virtual double c_star(){ return imp_->c_star(); }

  protected:
    ParamImp *imp_;
};

}

#endif
