#ifndef LANGEVIN_PATH_SCHEME_H
#define LANGEVIN_PATH_SCHEME_H

#include <iostream>

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
        void trajectory( gsl_rng *r, ParameterType &c, PathType &out );
        ComplexType sample_transition_density(gsl_rng *r, ComplexType c );
        void setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out );
        void generate_random_starts( gsl_rng *r, PathType &out );

        double log_prior_ratio(ParameterType &c_star, ParameterType &c);
        double log_transition_density_ratio(ParameterType &c_star, ParameterType &c);
        double log_path_likelihood( PathType &x, ParameterType &c, CoarsePathType &y );

        static ParameterType default_parameters( Options<LangevinDynamics>& o )
        {
            int M;
            M = o.cutoff();
            int D = 2*M*(M+1);
            ParameterType real_c(D);
            //real_c(D-1) = ComplexType(0.5, -0.5);
            return real_c;
        }

    
    protected:
        Options<LangevinDynamics> opts;
        int _cutoff;
        int _parameter_dimension;


}; // end of class LangevinDynamics

}
#endif
