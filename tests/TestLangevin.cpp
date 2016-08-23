#include "LangevinDynamics.h"

#define DEBUG
#include "Options.h"

#include "lest/lest.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



using namespace std;
using namespace MCMC;
using lest::approx;

const lest::test specification[] =
{

    CASE( "LangevinDynamics_log_transition" )
    {   
        Options opts;
        LangevinDynamics ld( opts );
        LangevinDynamics::ParameterType c_star(4);
        LangevinDynamics::ParameterType c(4);
        double result = ld.log_transition( c_star, c );
        EXPECT( 0 == result );
    },

    CASE( "LangevinDynamics_log_path_likelihood" )
    {
        Options opts;
        LangevinDynamics ld( opts );
        opts.set_observation_noise_sigma( 0.5 );
        
        LangevinDynamics::PathType x(1, 1, 2, 2);
        // Set up the real path.
        x(0,0,0,0) = 0.5454;
        x(0,0,0,1) = 0.8;
        x(0,0,1,0) = 0.5454;
        x(0,0,1,1) = 0.8;
        x(0,0,0,0) = 0.56;
        x(0,0,0,1) = 0.73;
        x(0,0,1,0) = 0.52;
        x(0,0,1,1) = 0.71;
        
        LangevinDynamics::CoarsePathType y(1,2,2);
        // Set up the observations.
        y(0, 0, 0) = 0.54;
        y(0, 0, 1) = 0.8;
        y(0, 1, 0) = 0.5;
        y(0, 1, 1) = 0.7;
        
        LangevinDynamics::ParameterType c(4);
        // Set up the modes.
        c(0) = c(1) = c(2) = 0.0;
        c(3) = ComplexType(0.5, -0.5);
        
        // Diffusion coefficient
        double log_sigma = log(0.01);
        double expected;
        double result;

        expected = -26.5;
        result = ld.log_path_likelihood( x, c, log_sigma, y, 0, 0, 0);
        EXPECT( expected == approx(result));
        
        expected = -7586.99;
        result = ld.log_path_likelihood( x, c, log_sigma, y, 0, 0, 1);
        EXPECT( expected == approx(result));
        
        expected = -2.5;
        result = ld.log_path_likelihood( x, c, log_sigma, y, 0, 1, 0);
        EXPECT( expected == approx(result));
    },
    
    CASE( "LangevinDynamics_forward_sim" )
    {
        Options opts;
        LangevinDynamics ld( opts );

        LangevinDynamics::CoarsePathType y(1,1,1);
        
        opts.set_observation_noise_sigma( 0.5 );
        double log_sigma = log(0.01);
        

        LangevinDynamics::ParameterType c(4);
        // Set up the modes.
        c(0) = c(1) = c(2) = 0.0;
        c(3) = ComplexType(0.5, -0.5);
               
        int K = 1;
        int L = 2;
        int M = 8;
        
        LangevinDynamics::PathType x(K, L, M, 2);
        opts.set_extra_data_ratio(M);
        
        EXPECT( opts.extra_data_ratio() == M );
        
        gsl_rng *r;
        gsl_rng_env_setup();
        r = gsl_rng_alloc( gsl_rng_default );
        gsl_rng_set( r, opts.rng_seed() );
        
        for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L-1; ++l )
        {
        for( size_t m=1; m<M; ++m )
        {
        ld.forward_sim(r, c, log_sigma, y, k, l, m, x);   
        }
        ld.forward_sim(r, c, log_sigma, y, k, l+1, 0, x);
        }

         for( size_t k=0; k<K; ++k )
         for( size_t l=0; l<L; ++l )
         {
         std::cout<<"+ "<< x(k,l,0,0) << " " << x(k,l,0,1) << std::endl; 
         if( l<L-1 )
         for( size_t m=1; m<M; ++m )
            std::cout<< ". " << x(k,l,m,0) << " " << x(k,l,m,1) << std::endl; 
         }
    },
    
    CASE( "LangevinDynamics_log_prior" )
    {
        Options opts;
        LangevinDynamics ld( opts );
        opts.set_parameter_proposal_sigma( 1.0 );        
        LangevinDynamics::ParameterType c(4);
        // Set up the modes.
        c(0) = c(1) = c(2) = 0.0;
        c(3) = ComplexType(0.5, -0.5);

        double result = ld.log_prior(c);
        double expected = -625.0;
        EXPECT(result == approx(expected) );
    },
    
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv , std::cout  );
}
