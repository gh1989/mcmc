#ifndef LANGEVIN_DYNAMICS_H
#define LANGEVIN_DYNAMICS_H

#include <string>
#include <fstream>
#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Dynamics.h"
#include "FourierSeries.h"
#include "Options.h"
using namespace MCMC;

namespace MCMC
{

//A basic one mode potential Langevin class
class LangevinDynamics : public DynamicsBase {
    
    public:    
    
        typedef ComplexType ParameterPointType;    
        typedef Tensor<ParameterPointType, 1> ParameterType;
        typedef Tensor<ParameterPointType, 2> ParameterChainType;
        
        typedef Tensor<double, 3> CoarsePathType;
        typedef Tensor<double, 4> PathType;    
        typedef Tensor<double, 1> SigmaChainType;
        
        LangevinDynamics(Options &o) : _opts(o), DynamicsBase(o)
        {
             std::cout<<  "LangevinDynamics(Options &o) called in " << this << std::endl;
             _cutoff = _opts.cutoff();
             _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }

        // Do not instantiate without options.
        LangevinDynamics() = delete;
        
        LangevinDynamics(const LangevinDynamics& other) : _opts(other.opts()), DynamicsBase(other.opts())
        {
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }
         
        LangevinDynamics(const LangevinDynamics&& other): _opts(other.opts()), DynamicsBase(other.opts())
        {
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }
 
        LangevinDynamics& operator=(const LangevinDynamics& other)
        {
            _opts = other.opts();
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
            return *this;
        }
        
        LangevinDynamics& operator=(const LangevinDynamics&& other)
        {
            _opts = other.opts();
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
            return *this;
        }
        
        ~LangevinDynamics()
        {}
        
        int parameter_dimension()
        {
            return _parameter_dimension; 
        }
        
        int parameter_dimension( const Options& o )
        {
            _cutoff = _opts.cutoff();
            return 2*_cutoff*(_cutoff+1);
        }

        static std::string dynamics_string()
        {
            return "Langevin_";
        }
        
        void output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file);
        ComplexType sample_transition_density(gsl_rng *r, ComplexType c );
        ParameterType sample_transition_density(gsl_rng *r, ParameterType& C);
        void forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out );
        double log_prior(ParameterType &c);
        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_path_likelihood( PathType &x, ParameterType &c, double sigma, CoarsePathType &y, int k, int l, int m );

        static ParameterType default_parameters( Options& o )
        {
            int M;
            M = o.cutoff();
            int D = 2*M*(M+1);
            ParameterType real_c(D);
            for(size_t i=0; i<D; ++i)
                real_c(i) = ComplexType(0, 0);
            
            // If we are not inferring the drift then turn off potential as a preference.
            bool have_a_potential = o.infer_drift_parameters();
            bool single_mode = o.single_mode();
            if (have_a_potential)
            {
                
                real_c(D-1) = ComplexType(0.5, -0.5);
            
                if ( (M >= 1) && (!single_mode) )
                {
                    real_c(0) = ComplexType(0.5,    -0.5);
                    real_c(1) = ComplexType(0.25,    0.5);
                    real_c(2) = ComplexType(1.0,     1.0);
                    real_c(3) = ComplexType(-1.0,   -1.7);
                }
            
                if ( (M >= 2) && (!single_mode) )
                {
                    real_c(4) = ComplexType(0.5,   -0.5);
                    real_c(5) = ComplexType(0.5,   -0.5);
                    real_c(6) = ComplexType(0.25,    1.0);
                    real_c(7) = ComplexType(0.15,    1.0);
                    real_c(8) = ComplexType(0.5,   -0.5);
                    real_c(9) = ComplexType(0.5,   -0.5);
                    real_c(10) = ComplexType(0.7,   0.2);
                    real_c(11) = ComplexType(0.2,   0.7);
                }
                
            }
            
            return real_c;
        }
        
    
    protected:
        Options _opts;
        int _cutoff;
        int _parameter_dimension;


}; // end of class LangevinDynamics

}
#endif
