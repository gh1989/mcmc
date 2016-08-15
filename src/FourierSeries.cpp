#include <iostream>
#include <iomanip>
#include "FourierSeries.h"

using namespace MCMC;

FourierSeries::FourierSeries( int M_ ) : M(M_), total_modes( (2*M_+1)*(2*M_+1) )
{
    modes = Tensor<ComplexType, 2>( 2*M+1, 2*M+1 );
    
    for(size_t i=0; i<2*M+1; ++i)
        for(size_t j=0; j<2*M+1; ++j)
            modes(i,j) = ComplexType(0,0);
}
     
void FourierSeries::set_mode( int i, int j, ComplexType mode)
{
    int ii = M + i;
    int jj = M + j;
    modes(ii , jj ) = mode;

    // The reality constraint.
    ii = M - i;
    jj = M - j;
    modes(ii , jj ) = std::conj(mode);
}

void FourierSeries::set_modes( Tensor<ComplexType, 1> &c )
{
    int idx;

    int ii;
    int jj; 
   
    for( size_t i=1; i<M+1; ++i )
        set_mode( i, 0, c(i-1) );

    for( int i=-M; i<M+1; ++i )
        for( size_t j=1; j<M+1; ++j )
        { 
            idx = M + (j-1)*(2*M+1) + (i+M);
            set_mode( i, j, c( idx ) );
        }
}

ComplexType FourierSeries::get_mode( int i, int j )
{
    if (i==0 && j==0) 
        return ComplexType(0,0);

    int ii = M + i;
    int jj = M + j;
    return modes(ii, jj);
}

Vector2d FourierSeries::grad( double x, double y )
{    
    Vector2d ret;
    
    ComplexType x_component(0,0);
    ComplexType y_component(0,0);    
    ComplexType mode;
    ComplexType tmp;
    
    for( int i= -M; i<M+1; ++i )
    for( int j= -M; j<M+1; ++j )
    {
        mode = get_mode( i, j );
        tmp = mode * std::exp( _2PI_I*(i*x+j*y) );    
        x_component += ComplexType(i,0)*_2PI_I*tmp;
        y_component += ComplexType(j,0)*_2PI_I*tmp;
    }

    ret(0) = std::real( x_component );
    ret(1) = std::real( y_component );
    
    return ret;
}

Vector2d FourierSeries::grad( Vector2d v )
{    
    return grad( v(0), v(1) );
}

double FourierSeries::evaluate( double x, double y )
{
    ComplexType ret(0, 0);
    ComplexType tmp;
    ComplexType mode;
    
    for( int i=-M; i<M+1; ++i )
        for( int j=-M; j<M+1; ++j )
        {
            mode = get_mode( i, j );
            ret += mode*std::exp( _2PI_I*(i*x+j*y) );
        }        
    
    return std::real( ret );
}

double FourierSeries::evaluate( Vector2d v )
{
    return evaluate( v(0), v(1) );
}

/*
 * Helper functions.
 */

void FourierSeries::print_modes()
{
    printf("Potential modes: \n");
    ComplexType mode;
    for( int j=M; j>-M-1; --j )
    {
        for( int i=-M; i<M+1; ++i )
        {
            mode = get_mode(i,j);
            std::cout<< std::fixed;
            std::cout<< std::setprecision(2);
            std::cout << mode << "\t\t";
        }
        std::cout<< std::endl;
    }
}
