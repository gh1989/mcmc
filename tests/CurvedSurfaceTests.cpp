#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>
#include <string>

#include "CurvedSurfaceDynamics.h"
#include "Options.h"


using namespace std;

int main( int argc, char *argv[] )
{
    Options o( argc, argv );
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, o.rng_seed() );

    CurvedSurfaceDynamics curved_dyn(o);
    
    Vector2d x;
    x << 1.0, 1.0;
    
    CurvedSurfaceDynamics::ParameterType c(4);
    
    // Try it out with a flat surface, the drift should be zero?
    for(size_t i=0; i<4; ++i)
        c(i) = 0.0;    
    
    std::cout << curved_dyn.drift( x, c ) << std::endl;
    std::cout << curved_dyn.diffusion( x, c ) << std::endl;
    
    // Now try with a simple h function.
    c(0) = 0.0;   
    c(1) = 0.0;   
    c(2) = 0.0;   
    c(3) = ComplexType(0.5,-0.5);   

    std::cout << curved_dyn.drift( x, c ) << std::endl;
    std::cout << curved_dyn.diffusion( x, c ) << std::endl;

    double delta_x = 0.001;
    size_t number_points = 2000;
 
    
    return 42;
}
