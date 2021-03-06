#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath> 

#include "BridgeDynamics.h"

using namespace MCMC;

namespace MCMC
{
    
template<class Dynamics_>
class Particle
{
public:

    using ParameterType         = typename Dynamics_::ParameterType;
    using ParameterChainType    = typename Dynamics_::ParameterChainType;
    using PathType              = typename Dynamics_::PathType;
    using CoarsePathType        = typename Dynamics_::CoarsePathType;

   Particle(Options& o) : _dynamics(o)
    {
        K = o.parallel_paths();
        L = o.path_length();
        M = o.extra_data_ratio();
        x = PathType( K, L, M, 2 ); 
    }
    
    // Do not instantiate without options.
    Particle() = delete;
    
    Particle( const Particle<Dynamics_>& other ) : _dynamics(other.dynamics())
    {
        x = other.path();
        K = other.parallel_paths();
        L = other.path_length();
        M = other.extra_data_ratio();
    }
    
    Particle<Dynamics_>& operator=(const Particle<Dynamics_>& other)
    {
        _dynamics = other.dynamics();
        x = other.path();
        K = other.parallel_paths();
        L = other.path_length();
        M = other.extra_data_ratio();
        return *this;
    }
    
    Particle(const Particle<Dynamics_>&& other) : _dynamics(other.dynamics())
    {
        x = other.path();
        K = other.parallel_paths();
        L = other.path_length();
        M = other.extra_data_ratio();
    }
    
    Particle<Dynamics_>& operator=(const Particle<Dynamics_>&& other)
    {
        _dynamics = other.dynamics();
        x = other.path();
        K = other.parallel_paths();
        L = other.path_length();
        M = other.extra_data_ratio();
        return *this;
    }
    
    ~Particle(){}

    double unnormal_weight( size_t t, ParameterType &c, double sigma, CoarsePathType &y);
    void setup_starts(gsl_rng *r, CoarsePathType &y);
    void forward_sim(gsl_rng *r, ParameterType &c, double sigma, CoarsePathType &y, size_t t );
    const PathType& path() const { return x; }
    const Dynamics_& dynamics() const { return _dynamics; }
    const size_t parallel_paths() const { return K; }
    const size_t path_length() const  { return L; }
    const size_t extra_data_ratio() const { return M; }
    
private:  
    size_t K;
    size_t L;
    size_t M;
    PathType x;
    Dynamics_ _dynamics;
};


template<class Dynamics_>
double Particle<Dynamics_>::unnormal_weight(size_t t, ParameterType &c, double sigma, CoarsePathType &y)
{   
    double log_total = 0;
    for( size_t k=0; k<K; ++k ) // p( y | x, c )
        log_total += _dynamics.log_p( y, x, c, sigma, k, t );
    return log_total;
}

template<>
double Particle<BridgeDynamics>::unnormal_weight(size_t t, ParameterType &c, double log_sigma, CoarsePathType &y)
{   
    double log_total = 0;
    double temp = 0; 
    for( size_t k=0; k<K; ++k )
    {
        /*
        for(size_t m=1; m<M+1; ++m)
        {
            // \hat{p}( x(n + dt) ... x(n+1) | y(n+1), x(n), c )
            temp = _dynamics.log_path_likelihood( x, c, log_sigma, y, k, t-1, m );
            std::cout << "numerator: " << temp << std::endl;
            log_total += temp;
            
            // \hat{p}( x(n + dt) ... x(n+1) | y(n+1), x(n), c )
            temp = _dynamics.log_bridge_likelihood( x, c, log_sigma, y, k, t-1, m ); 
            std::cout << "denominator: " <<temp << std::endl;
            log_total -= temp;
        } 
        */
        // p( y | x, c )
        // log_total += _dynamics.log_path_likelihood( x, c, log_sigma, y, k, t, 0 );
        
        // p( y | x, c )
        temp= _dynamics.log_p( y, x, c, log_sigma, k, t );        
        if (std::isnan(temp))
        {
            std::cout<< "WARNING! temp is nan?!" << std::endl;
            if(k!=0)
                std::cout<< "k is not zero: " << k << std::endl;       
        }
        log_sigma += temp;
    }        

    std::cout<<"Particle<BridgeDynamics>::unnormal_weight log_total:" << log_total << std::endl;
    
    if (std::isnan(log_total))
    {
        std::cout<< "WARNING! log_total is nan?!" << std::endl;
        std::cout<< "t" << t << std::endl;
        std::cout<< "y" << y << std::endl;
        std::cout<< "c" << c << std::endl;
        std::cout<< "log_sigma" << log_sigma << std::endl;        
    }
    
    return log_total;
}

template<class Dynamics_>
void Particle<Dynamics_>::setup_starts(gsl_rng *r, CoarsePathType &y)
{
    for( size_t k=0; k<K; ++k )
    {
        x(k, 0, 0, 0 ) = y(k, 0, 0);
        x(k, 0, 0, 1 ) = y(k, 0, 1);
    }
}

template<>
void Particle<BridgeDynamics>::forward_sim(gsl_rng *r, ParameterType &c, double sigma, CoarsePathType &y, size_t l )
{
    /*
    Forward simulate from (l-1,l]
    */    
    for( size_t k=0; k<K; ++k )
    {
        for( size_t m=1; m<M; ++m )
            _dynamics.forward_sim(r, c, sigma, y, k, l-1, m, x);   
        _dynamics.forward_sim(r, c, sigma, y, k, l-1, M, x);
    }  
}

template<class Dynamics_>
void Particle<Dynamics_>::forward_sim(gsl_rng *r, ParameterType &c, double sigma, CoarsePathType &y, size_t l )
{
    /*
    Forward simulate from (l-1,l]
    */    
    for( size_t k=0; k<K; ++k )
    {
        for( size_t m=1; m<M; ++m )
            _dynamics.forward_sim(r, c, sigma, y, k, l-1, m, x);   
        _dynamics.forward_sim(r, c, sigma, y, k, l-1, M, x);
    }  
}

}
#endif
