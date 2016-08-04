#ifndef OU_DYNAMICS_H
#define OU_DYNAMICS_H

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "Dynamics.h"

using namespace MCMC;

namespace MCMC
{

class OUDynamics : public DynamicsBase
{
    public:
        typedef Tensor<double, 1> ParameterType;
        typedef Tensor<double, 2> ParameterChainType;

        typedef Tensor<double, 3> CoarsePathType;
        typedef Tensor<double, 4> PathType;    

        OUDynamics(Options<OUDynamics>& o): opts(o)
        {
            _d_sigma = opts.diffusion_coefficient();  
            _dt = opts.trajectory_path_delta();
            _d_variance = _d_sigma*_d_sigma*_dt;
            _d_const = (0.5 / _d_variance);
            _o_variance = opts.observation_noise_variance();
            _o_const = (0.5 / _o_variance);
            _parameter_dimension = 1;
        }

        OUDynamics()
        {
            _parameter_dimension = 1;
        }

        int parameter_dimension(){ return _parameter_dimension; }
        int parameter_dimension( const Options<OUDynamics>& o ){ return _parameter_dimension; }

        void output_file_timeseries(ParameterChainType &ccc);
        double sample_transition_density(gsl_rng *r, double c );

        void forward_sim( gsl_rng *r, ParameterType &c,
                          int k, int l, int m, 
                          PathType &out );        

        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_prior(ParameterType &c);
        double log_path_likelihood( PathType &x, 
                                    ParameterType &c, 
                                    CoarsePathType &y,
                                    int k, int l, int m );

        static ParameterType default_parameters( Options<OUDynamics>& o )
        {
            ParameterType real_c(1);
            real_c(0) = 1.70;
            return real_c;
        }


    protected:
        Options<OUDynamics> opts;

        double _d_sigma;
        double _dt;
        double _d_variance; 
        double _d_const;
        double _o_variance;
        double _o_const;
        
        int _parameter_dimension;

}; // end of class OUDynamics

}
#endif

