#ifndef FOURIER_SERIES_H
#define FOURIER_SERIES_H

#include <complex>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

typedef std::complex<double> ComplexType;

const ComplexType _2PI_I = 2*M_PI*ComplexType(0,1);

namespace MCMC
{

class FourierSeries
{
public:
    FourierSeries(int M_);
    //FourierSeries(){};
    ~FourierSeries();

    Vector2d grad( Vector2d );
    Vector2d grad( double x, double y );

    double evaluate( Vector2d );
    double evaluate( double x, double y );

    void set_mode( int i, int j, ComplexType );
    ComplexType get_mode( int, int );

    void print_modes();

private:
    int M;
    int total_modes;
    Tensor<ComplexType, 2> modes;
};

}

#endif
