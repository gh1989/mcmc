#include "FourierSeries.h"
#include "lest/lest.hpp"

using namespace MCMC;

const lest::test specification[] =
{

    CASE( "FourierSeries_operator")
    {
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
    
        Tensor<ComplexType, 1> c(4);
        int cutoff = 1;
        FourierSeries f(cutoff);
        int idx;
        int ii;
        int jj; 
        
        for( size_t i=1; i<cutoff+1; ++i )
            f.set_mode( i, 0, c(i-1) );

        for( size_t i=-cutoff; i<cutoff+1; ++i )
            for( size_t j=1; j<cutoff+1; ++j )
            {
                ii  = cutoff + i; 
                jj  = cutoff + j;
                idx = ii + jj*(2*cutoff+1);
                f.set_mode( i, j, c( idx ) );
            }

        EXPECT( f.get_mode(  0, 0 ) == 0.0  );
        EXPECT( f.get_mode(  1, 0 ) == c(0) );
        EXPECT( f.get_mode( -1, 1 ) == c(1) );
        EXPECT( f.get_mode(  0, 1 ) == c(2) );
        EXPECT( f.get_mode(  1, 1 ) == c(3) );
    },

};


int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv /*, std::cout*/  );
}
