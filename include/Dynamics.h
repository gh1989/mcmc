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
    DynamicsBase(Options &o)
    {
        _dt = opts.trajectory_path_delta();
        _o_variance = opts.observation_noise_variance();
        _o_const = (0.5 / _o_variance);
    }
    DynamicsBase(){}
    ~DynamicsBase(){}
    
protected:
    double _dt;
    double _o_variance;
    double _o_const;
    Options opts;
    int _parameter_dimension;
};


}
#endif // DYNAMICS_H
