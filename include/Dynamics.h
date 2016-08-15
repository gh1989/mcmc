#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"

namespace MCMC
{

class DynamicsBase
{
public:
    DynamicsBase(Options &o) : _opts(o) {}
    DynamicsBase(){}
    ~DynamicsBase(){}
    
    virtual double log_prior_sigma(double log_sigma);
    virtual double log_transition_sigma(double log_sigma_star, double log_sigma);
    
    const Options &opts() const {return _opts;}
    
protected:
    Options _opts;
    int _parameter_dimension;
};

}

#endif // DYNAMICS_H
