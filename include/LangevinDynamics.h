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

#include "Dynamics.h"
#include "FourierSeries.h"
#include "Options.h"
using namespace MCMC;

namespace MCMC
{

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
        { 
            //std::cout<<"~LangevinDynamics() called in "<< this << std::endl;
        }
        
        int parameter_dimension()
        {
            return _parameter_dimension; 
        }
        
        int parameter_dimension( const Options& o )
        {
            _cutoff = _opts.cutoff();
            return 2*_cutoff*(_cutoff+1);
        }

        void output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss);
        ComplexType sample_transition_density(gsl_rng *r, ComplexType c );
        void forward_sim( gsl_rng *r, ParameterType &c, double sigma, int k, int l, int m, PathType &out );
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
            real_c(D-1) = ComplexType(0.5, -0.5);
            return real_c;
        }
        
    
    protected:
        Options _opts;
        int _cutoff;
        int _parameter_dimension;


}; // end of class LangevinDynamics

}
#endif
