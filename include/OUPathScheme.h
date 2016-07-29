#ifndef OU_SCHEME_H
#define OU_SCHEME_H

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Options.h"
#include "OUParams.h"
#include "PathScheme.h"

using namespace MCMC;

namespace MCMC
{

class OUPathSchemeImp : public PathSchemeImp
{
    public:
        OUPathSchemeImp(Options &opts_): PathSchemeImp(opts_), opts(opts_)
        {
            _params = OUParams(opts_);
  
            int L = opts.path_length();
            int K = opts.parallel_paths();
            int M = opts.extra_data_ratio();

            current     = Tensor<double,4>( K, L, M, 2 );
            proposed    = Tensor<double,4>( K, L, M, 2 );
            observed    = Tensor<double,3>( K, L, 2 );
            real        = Tensor<double,4>( K, L, M, 2 );
        }
        
    /* virtual */
    void    generate_true_trajectory( gsl_rng *r );
    void    generate_observations( gsl_rng *r );
    void    propose_trajectory( gsl_rng *r );
    void    propose_parameters( gsl_rng *r ) { params().propose_parameters(r); }
    void    accept() { params().accept(); }
    void    finish(){ params().output_file_timeseries(); }
    void    store_chain(int n) { params().store_chain(n); }
    double  log_path_ratio();
    double  log_prior_ratio(){ return params().log_prior_ratio(); }
    double  log_transition_density_ratio(){ return params().log_transition_density_ratio(); } 
    /* end virtual */

    //////////////////////////////////////////////////////////////////////////

    void ornstein_uhlenbeck_trajectory( gsl_rng *r, double c, Tensor<double, 4> &out );
    void setup_observed_starts( gsl_rng *r, const Tensor<double, 3> &y, Tensor<double, 4> &out );
    void generate_random_starts( gsl_rng *r, Tensor<double, 4> &out );

    OUParams params() { return _params; }

    protected:
        Options opts;
        OUParams _params;

        static constexpr double true_c = 1.70; // Hardcoded

        Tensor<double, 4> current;
        Tensor<double, 4> proposed;
        Tensor<double, 3> observed;
        Tensor<double, 4> real;       

}; // end of class OUPathSchemeImp

class OUPathScheme : public PathScheme
{
    public:
        OUPathScheme(Options &opts_) {
            imp_ = new OUPathSchemeImp(opts_);
        }
}; // end of class OUPathScheme

}
#endif

