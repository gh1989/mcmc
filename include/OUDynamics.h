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
        typedef double ParameterPointType;    
        typedef Tensor<ParameterPointType, 1> ParameterType;
        typedef Tensor<ParameterPointType, 2> ParameterChainType;

        typedef Tensor<double, 3> CoarsePathType;
        typedef Tensor<double, 4> PathType;    
        
        OUDynamics(Options& o): opts(o), DynamicsBase(o)
        {
            _parameter_dimension = 1;
        }

        OUDynamics()
        {
            _parameter_dimension = 1;
        }

        int parameter_dimension(){ return _parameter_dimension; }
        int parameter_dimension( const Options& o ){ return _parameter_dimension; }

        void output_file_timeseries(ParameterChainType &ccc);
        double sample_transition_density(gsl_rng *r, double c );

        void forward_sim( gsl_rng *r, 
                          ParameterType &c,
                          double sigma,
                          int k, int l, int m, 
                          PathType &out );    

        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_prior(ParameterType &c);
        double log_path_likelihood( PathType &x, 
                                    ParameterType &c, 
                                    double sigma,
                                    CoarsePathType &y,
                                    int k, int l, int m );

        static ParameterType default_parameters( Options& o )
        {
            ParameterType real_c(1);
            real_c(0) = 1.70;
            return real_c;
        }


    protected:
        Options opts;

}; // end of class OUDynamics

}
#endif

