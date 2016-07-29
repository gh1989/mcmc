#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LikelihoodFreeScheme.h"
#include "Options.h"

namespace MCMC
{

class Algorithm
{
    public:
        Algorithm( const Options& opts, LikelihoodFreeScheme& algo_scheme_ ) : algo_scheme(algo_scheme_), _opts(opts) {};
        ~Algorithm(){};
        void run( gsl_rng *r );

        Options opts(){ return _opts; }

    private:
        Options _opts;
        LikelihoodFreeScheme algo_scheme;
};

}
#endif
