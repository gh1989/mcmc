#include "FourierSeries.h"
#include "lest/lest.hpp"

using namespace MCMC;
using lest::approx;

const lest::test specification[] =
{

    CASE( "FourierSeries_operator")
    {
        approx eps = approx::custom().epsilon( 1e-100 );
        
        FourierSeries f(1);
        FourierSeries f_new(1);
        f.set_mode(1,1,0.25);
        f_new = f;
        EXPECT( f_new.get_mode(1,1) == 0.25 );

        FourierSeries f_newer;
        f_newer = f;
        EXPECT( f_newer.get_mode(1,1) == 0.25 );
    },

    CASE( "FourierSeries_set_mode")
    {               
        int cutoff = 1;
        FourierSeries f(cutoff);
        
        f.set_mode(0, 1, 1.0);
        EXPECT(f.get_mode(0,1) == 1.0 );
        Tensor<ComplexType, 1> c(4);
        for(size_t i=0; i<4; ++i) 
            c(i) = i;
        std::cout << "Parameters";
        std::cout << c << std::endl;
        f.set_modes(c);                   
                   
        EXPECT( f.get_mode(  0, 0 ) == 0.0 );
        EXPECT( f.get_mode(  1, 0 ) == 0.0 );
        EXPECT( f.get_mode(  -1, 1 ) == 1.0 );
        EXPECT( f.get_mode(  0, 1 ) == 2.0 );
        EXPECT( f.get_mode(  1, 1 ) == 3.0 );

        f.set_mode( 1,1, ComplexType(0,1) );
        EXPECT( f.get_mode(-1,-1) == ComplexType(0,-1) );
        
    },
    
    CASE("FourierSeries_grad")
    {
        Vector2d grad_result;
        
        int cutoff = 1;
        FourierSeries f(cutoff);
        
        f.set_mode(0, 1, 1.0);
        EXPECT(f.get_mode(0,1) == 1.0 );
        Tensor<ComplexType, 1> c(4);
        for(size_t i=0; i<3; ++i) 
            c(i) = 0;
        c(3) = ComplexType(0.5, -0.5);
        f.set_modes(c);                 
        grad_result = f.grad( 0.5454, 0.8);
        
        std::cout<< grad_result;
        
        EXPECT( -8.73253 == approx(grad_result(0)) );
        EXPECT( -8.73253 == approx(grad_result(1)) );
        
        Vector2d x(0.5454, 0.8);
        
        EXPECT( f.grad(x) == grad_result );
        
        }

};


int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv /*, std::cout*/  );
}
