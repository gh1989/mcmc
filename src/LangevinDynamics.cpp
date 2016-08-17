#include "LangevinDynamics.h"

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

ComplexType LangevinDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    double random_real_part = gsl_ran_gaussian(r, proposal_sigma);
    double random_imag_part = gsl_ran_gaussian(r, proposal_sigma);
    return ComplexType(random_real_part, random_imag_part) + c;
}

double LangevinDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    return 0;
}

double LangevinDynamics::log_prior(ParameterType &c)
{
    double log_total = 0;
    double parameter_sigma = _opts.parameter_proposal_sigma();
    double exponential_constant = 0.5 / (parameter_sigma * parameter_sigma);
    int c_dim = parameter_dimension();

    for( size_t i=0; i<c_dim; ++i)
        log_total -= exponential_constant*( pow( std::abs(c(i)), 2) );

    return log_total;
}

double LangevinDynamics::log_path_likelihood(  PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m )
{
    double log_total = 0;
    double dt = _opts.trajectory_path_delta();
    double obs_sigma = _opts.observation_noise_sigma();
    FourierSeries V(_cutoff);
    V.set_modes( c );

    Vector2d xt;
    Vector2d xtminus1;
    Vector2d gradvtminus1;
    
    // dX_t = -grad( V ) dt + sqrt(2D) dW_t
    // D = Exp( log_sigma )
    
    double exponential_constant_path = 0.25 / ( exp(2*log_sigma)*dt );
    double exponential_constant_obs = 0.5 / ( obs_sigma*obs_sigma );
    
    if( m == 0 )
    {        
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    {
        xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
        gradvtminus1 = V.grad( xtminus1 );
        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2);
        log_total -= exponential_constant_path * pow( xt(1)-xtminus1(1)+gradvtminus1(1)*dt, 2);
    }
    
    return log_total;
}

void LangevinDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, int k, int l, int m, PathType &out )
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);
    
    FourierSeries V(_cutoff);
    V.set_modes(c);
            
    Vector2d xtminus1;
    Vector2d gradtminus1;
    
    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);

    gradtminus1 = V.grad( xtminus1 );
 
    out(k, l, m, 0 ) = xtminus1(0)-gradtminus1(0)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));
    out(k, l, m, 1 ) = xtminus1(1)-gradtminus1(1)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));  
}

void LangevinDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss)
{
    cout<<"LangevinDynamics::output_file_timeseries"<<endl;
    size_t N  = _opts.mcmc_trials();
    
    unsigned long milliseconds_since_epoch = 
    std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
        
    std::string filename = "output/";
    filename += _opts.output_subfolder();
    std::cout << "Putting in subfolder:" << _opts.output_subfolder() << std::endl;
    filename += "/LangevinTimeSeries_";
    filename += std::to_string( milliseconds_since_epoch );
    filename += ".txt";

    std::cout<< "Output file: " << filename << std::endl;
    
    std::cout<< "Currently in " << this << "- LangevinDynamics::output_file_timeseries" << std::endl;
    std::cout<< "Options are: " << std::endl;
    _opts.print_options(std::cout);
    
    ofstream mcmc_file;
    mcmc_file.open(filename);

    int defining_modes = parameter_dimension(_opts);

    // experiment options
    _opts.print_header( mcmc_file );
    
    // headers
    mcmc_file << "#";
    for( size_t m=0; m<defining_modes; ++m)             
        mcmc_file << m << "\t" << m << "\t";
    mcmc_file << "sigma" << std::endl;
    
    // data
    for( size_t n=0; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(ccc(n,m)) << "\t" << std::imag(ccc(n,m)) << "\t";
        mcmc_file << sss(n) << std::endl;
    }
    mcmc_file.close();
} 
