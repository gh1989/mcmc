#include "OUDynamics.h"
#include "Options.h"
#include "lest/lest.hpp"

using namespace std;
using namespace MCMC;

const lest::test specification[] =
{

   CASE( "OUDynamics_log_transition" )
    {   
        Options opts;
        OUDynamics ou( opts );
        OUDynamics::ParameterType c_star(4);
        OUDynamics::ParameterType c(4);
        double result = ou.log_transition( c_star, c );
        EXPECT( 0 == result );
    },

    CASE( "OUDynamics_log_path_likelihood" )
    {
        Options opts;
        OUDynamics ou( opts );
        opts.set_observation_noise_sigma( 0.5 );
        
        OUDynamics::PathType x(1, 1, 2, 2);
        // Set up the real path.
        x(0,0,0,0) = 1.0;
        x(0,0,0,1) = 1.1;
        x(0,0,1,0) = 1.0;
        x(0,0,1,1) = 1.1;
        
        OUDynamics::CoarsePathType y(1,1,2);
        // Set up the observations.
        y(0, 0, 0) = 1.0;
        y(0, 0, 1) = 0.0;

        OUDynamics::ParameterType c(1);
        // Set up the modes.
        c(0) = 1.0;
        
        // Diffusion coefficient
        log_sigma = log(0.01);

        double result = ou.log_path_likelihood( x, c, log_sigma, y, 0, 0, 0);
        
    }

};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv , std::cout  );
}
