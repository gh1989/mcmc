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
        typedef Tensor<double, 1> SigmaChainType;        
        
        OUDynamics(Options& o): _opts(o), DynamicsBase(o)  {}

        OUDynamics() = default;
        OUDynamics& operator=(OUDynamics& other) = delete;
        OUDynamics( OUDynamics&& ) = delete;
        OUDynamics( OUDynamics& ) = delete;
        OUDynamics& operator=(OUDynamics&& other)
        {
            _opts = other.opts();
        }
        
        
        
        int parameter_dimension(){ return _parameter_dimension; }
        int parameter_dimension( const Options& o ){ return _parameter_dimension; }

        void output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file);
        void forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out );    

        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_prior(ParameterType &c);
        double log_path_likelihood( PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m );
        double sample_transition_density(gsl_rng *r, double c );
        
        static std::string dynamics_string()
        {
            return "OU_";
        }
        
        static ParameterType default_parameters( Options& o )
        {
            ParameterType real_c(1);
            real_c(0) = 1.70;
            return real_c;
        }


    protected:
        Options _opts;
        constexpr static double _parameter_dimension = 1;
        
}; // end of class OUDynamics

}
#endif

