#ifndef LANGEVIN_DYNAMICS_H
#define LANGEVIN_DYNAMICS_H

#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "FourierSeries.h"
#include "Options.h"
#include "Dynamics.h"

using namespace MCMC;

namespace MCMC
{

class LangevinDynamics  : public DynamicsBase {
    public:
    
        typedef Tensor<ComplexType, 1> ParameterType;
        typedef Tensor<ComplexType, 2> ParameterChainType;

        typedef Tensor<double, 3> CoarsePathType;
        typedef Tensor<double, 4> PathType;    

        LangevinDynamics(Options<LangevinDynamics> &o): opts(o)
        {
             _d_sigma = opts.diffusion_coefficient();  
             _dt = opts.trajectory_path_delta();
             _d_variance = _d_sigma*_d_sigma*_dt;
             _d_const = (0.5 / _d_variance);
             _o_variance = opts.observation_noise_variance();
             _o_const = (0.5 / _o_variance);
             _cutoff = opts.cutoff();
             _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }

        LangevinDynamics(){}

        int parameter_dimension(){ return _parameter_dimension; }
        int parameter_dimension( const Options<LangevinDynamics>& o )
        {
            _cutoff = opts.cutoff();
            return 2*_cutoff*(_cutoff+1);
        }

        void output_file_timeseries(ParameterChainType &ccc);
        ComplexType sample_transition_density(gsl_rng *r, ComplexType c );
        void forward_sim( gsl_rng *r, ParameterType &c, int k, int l, int m, PathType &out );
        double log_prior(ParameterType &c);
        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_path_likelihood( PathType &x, 
                                    ParameterType &c, 
                                    CoarsePathType &y,
                                    int k, int l, int m );

        static ParameterType default_parameters( Options<LangevinDynamics>& o )
        {
            int M;
            M = o.cutoff();
            int D = 2*M*(M+1);
            ParameterType real_c(D);
            for(size_t i=0; i<D; ++i)
                real_c(i) = ComplexType(0, 0);
            real_c(D-1) = ComplexType(0.5, -0.5);
            return real_c;
        }

    
    protected:
        Options<LangevinDynamics> opts;
        int _cutoff;
        int _parameter_dimension;
        double _d_sigma;
        double _dt;
        double _d_variance; 
        double _d_const;
        double _o_variance;
        double _o_const;


}; // end of class LangevinDynamics

}
#endif
