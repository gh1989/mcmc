#ifndef LANGEVIN_PATH_SCHEME_H
#define LANGEVIN_PATH_SCHEME_H

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include "Options.h"
#include "PathScheme.h"
#include "LangevinParams.h"

using namespace MCMC;

namespace MCMC
{

class LangevinPathSchemeImp : public PathSchemeImp {
    public:
        LangevinPathSchemeImp(Options &opts_): opts(opts_), PathSchemeImp(opts_)
        {
            _params = LangevinParams(opts_);     
        
            int L = opts.path_length();
            int K = opts.parallel_paths();
            int M = opts.extra_data_ratio();

            current     = Tensor<double,4>( K, L, M, 2 );
            proposed    = Tensor<double,4>( K, L, M, 2 );
            observed    = Tensor<double,3>( K, L, 2 );
            real        = Tensor<double,4>( K, L, M, 2 );

            // Hard coded
            int cutoff = 1;
            true_potential = FourierSeries(cutoff);
            true_potential.set_mode( 1, 1, ComplexType(0.5,-0.5) );
    
        }
        
    /* virtual */
    void    generate_true_trajectory( gsl_rng *r );
    void    generate_observations( gsl_rng *r );
    void    propose_trajectory( gsl_rng *r );
    void    propose_parameters( gsl_rng *r ) { params().propose_parameters(r); }
    void    accept() { params().accept(); }
    void    store_chain(int n) { params().store_chain(n); }
    double  log_path_ratio();
    double  log_prior_ratio(){ return params().log_prior_ratio(); }
    double  log_transition_density_ratio(){ return params().log_transition_density_ratio(); } 
    /* end virtual */

    //////////////////////////////////////////////////////////////////////////

    void langevin_trajectory( gsl_rng *r, FourierSeries &V, Tensor<double, 4> &out );
    void setup_observed_starts( gsl_rng *r, const Tensor<double, 3> &y, Tensor<double, 4> &out );
    void generate_random_starts( gsl_rng *r, Tensor<double, 4> &out );

    LangevinParams& params() { return _params; }

    protected:
        Options opts;
        LangevinParams _params;

        FourierSeries true_potential;

        Tensor<double, 4> current;
        Tensor<double, 4> proposed;
        Tensor<double, 3> observed;
        Tensor<double, 4> real;       

}; // end of class OUPathSchemeImp

class LangevinPathScheme : public PathScheme
{
    public:
        LangevinPathScheme(Options &opts_) {
            imp_ = new LangevinPathSchemeImp(opts_);
        }
}; // end of class OUPathScheme

}
#endif
