#ifndef BRIDGE_DYNAMICS_H
#define BRIDGE_DYNAMICS_H

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

class BridgeDynamics : public DynamicsBase {
    
    public:    
    
        typedef ComplexType ParameterPointType;    
        typedef Tensor<ParameterPointType, 1> ParameterType;
        typedef Tensor<ParameterPointType, 2> ParameterChainType;
        
        typedef Tensor<double, 3> CoarsePathType;
        typedef Tensor<double, 4> PathType;    
        typedef Tensor<double, 1> SigmaChainType;
        
        BridgeDynamics(Options &o) : _opts(o), DynamicsBase(o)
        {
             std::cout<<  "BridgeDynamics(Options &o) called in " << this << std::endl;
             _cutoff = _opts.cutoff();
             _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }

        // Do not instantiate without options.
        BridgeDynamics() = delete;
        
        BridgeDynamics(const BridgeDynamics& other) : _opts(other.opts()), DynamicsBase(other.opts())
        {
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }
         
        BridgeDynamics(const BridgeDynamics&& other): _opts(other.opts()), DynamicsBase(other.opts())
        {
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
        }
 
        BridgeDynamics& operator=(const BridgeDynamics& other)
        {
            _opts = other.opts();
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
            return *this;
        }
        
        BridgeDynamics& operator=(const BridgeDynamics&& other)
        {
            _opts = other.opts();
            _cutoff = _opts.cutoff();
            _parameter_dimension = 2*_cutoff*(_cutoff+1);
            return *this;
        }
        
        ~BridgeDynamics()
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
            return "Bridge_";
        }
        
        void output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file);
        ComplexType sample_transition_density(gsl_rng *r, ComplexType c );
        void forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out );
        double log_prior(ParameterType &c);
        double log_transition(ParameterType &c_star, ParameterType &c);
        double log_path_likelihood( PathType &x, ParameterType &c, double sigma, CoarsePathType &y, int k, int l, int m );
        double log_bridge_likelihood( PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m);

        static ParameterType default_parameters( Options& o )
        {
            int M;
            M = o.cutoff();
            int D = 2*M*(M+1);
            ParameterType real_c(D);
            for(size_t i=0; i<D; ++i)
                real_c(i) = ComplexType(0, 0);
            
            if (o.infer_drift_parameters())
                real_c(D-1) = ComplexType(0.5, -0.5);
            return real_c;
        }
        
    
    protected:
        Options _opts;
        int _cutoff;
        int _parameter_dimension;


}; // end of class BridgeDynamics

}
#endif
