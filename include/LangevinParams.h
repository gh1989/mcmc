#ifndef LANGEVIN_PARAMETERS_H
#define LANGEVIN_PARAMETERS_H

#include <Eigen/CXX11/Tensor>

#include "FourierSeries.h"

#include "Options.h"
#include "Params.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace MCMC;
using namespace Eigen;

namespace MCMC
{

class LangevinParamImp: public ParamImp {
    public:
        LangevinParamImp(Options &opts_) : _potential(FourierSeries(1)), _potential_star(FourierSeries(1)), ParamImp(opts_)
        {
            opts = opts_;
        
            // Hard coded for now.
            cutoff = 1;
            defining_modes = 2*cutoff*(cutoff+1);

            int N = opts.mcmc_trials();

            constants_chain = Tensor<ComplexType, 2>(N, defining_modes );
            diffusion_chain = Tensor<double, 3>(N, 2, 2);

            c       = Tensor<ComplexType, 1>( defining_modes );
            c_star  = Tensor<ComplexType, 1>( defining_modes );

            diffusion       = Matrix2d(2,2);
            diffusion_star  = Matrix2d(2,2);

            _potential      = FourierSeries(cutoff);
            _potential_star = FourierSeries(cutoff);

        }
        
        /* virtual */
        void propose_parameters(gsl_rng *r);
        double log_transition_density_ratio();
        double log_prior_ratio();
        void accept();
        void store_chain(int n);
        

        FourierSeries potential(){return _potential; }
        FourierSeries potential_star(){return _potential_star;}    
        /* end virtual */

    protected:
        Options opts;

        FourierSeries _potential;
        FourierSeries _potential_star;

        Matrix2d diffusion;
        Matrix2d diffusion_star;

        Tensor<ComplexType, 1> c;
        Tensor<ComplexType, 1> c_star;

        Tensor<ComplexType, 2> constants_chain;
        Tensor<double, 3> diffusion_chain;

        int cutoff;
        int defining_modes;
}; // end of class LangevinParamsImp

class LangevinParams: public Params {
  public:
    LangevinParams(Options &opts_) {
        imp_ = new LangevinParamImp(opts_);
    }
    LangevinParams(){}    
    
}; // end of class LangevinParams

}

#endif
