#ifndef OU_PARAMS_H
#define OU_PARAMS_H

#include <Eigen/CXX11/Tensor>

#include "Options.h"
#include "Params.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>

using namespace MCMC;
using namespace Eigen;

namespace MCMC
{

class OUParamImp: public ParamImp {
  public:
    OUParamImp(Options &opts_): ParamImp(opts_)
    {
        opts = opts_;
        int N = opts.mcmc_trials();
        std::cout<< "MCMC trials for OUParamImp: " << N << std::endl;
        constants_chain = Tensor<double, 1>(N);
        diffusion_chain = Tensor<double, 1>(N);

        _c = 0.0;
        _c_star = 0.0;

    }

    /* virtual */
    void propose_parameters(gsl_rng *r);
    double log_transition_density_ratio();
    double log_prior_ratio();
    void accept();
    void store_chain(int n);
    void output_file_timeseries();
    /* end virtual */
    
    double c(){ return _c; } 
    double c_star(){ return _c_star; }

  protected:
    Options opts;
    double _c;
    double _c_star;
    double sigma;
    double sigma_star;
    Tensor<double, 1> constants_chain;
    Tensor<double, 1> diffusion_chain;
};

class OUParams: public Params {
  public:
    OUParams(Options &opts_) {
        imp_ = new OUParamImp(opts_);
    }
    OUParams(){}
};

}

#endif
