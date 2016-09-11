#include "LangevinDynamics.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double LangevinDynamics::log_p( CoarsePathType &y,  PathType &x, ParameterType &c, double log_sigma, int k, int l)
{
    
    double obs_sigma =  _opts.observation_noise_sigma();
    double exponential_constant_obs  = 0.5 / ( pow( 2*M_PI*obs_sigma*obs_sigma, 1.5 ) );
   
    double log_total = 0;

    log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
    log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    
    return log_total;
}


ComplexType LangevinDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    double random_real_part = gsl_ran_gaussian(r, proposal_sigma);
    double random_imag_part = gsl_ran_gaussian(r, proposal_sigma);
    return ComplexType(random_real_part, random_imag_part) + c;
}

LangevinDynamics::ParameterType LangevinDynamics::sample_transition_density(gsl_rng *r, ParameterType& C)
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    bool single_mode = _opts.single_mode();
    size_t c_dim = parameter_dimension(_opts);
    ParameterType C_star(c_dim);
    for(size_t i=0; i<c_dim; ++i) C_star(i) = 0.0;
    int M = _opts.cutoff();

    int idx;
    
    double eigenvalue_power;
    if ((M==1)&&(single_mode))
    {
        eigenvalue_power = 2*M_PI*sqrt(2);
        C_star(c_dim-1) = ComplexType(  gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ), 
                                        gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ) ) + C(c_dim-1); 
    }
    else
    {
        for( size_t i=1; i<M+1; ++i )
        {
            eigenvalue_power = 2*M_PI*(i);
            C_star(i-1) = ComplexType( gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ), 
                                       gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ) ) + C(i-1); 
        }

        for( int i=-M; i<M+1; ++i )
            for( size_t j=1; j<M+1; ++j )
            { 
                idx = M + (j-1)*(2*M+1) + (i+M);
                eigenvalue_power = 2*M_PI*sqrt(i*i + j*j);
                C_star(idx) = ComplexType( gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ), 
                                           gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -1.0) ) ) + C(idx); 
            }
    }
    
    return C_star;
}

double LangevinDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    return 0;
}

// Gaussian prior: TODO: create a prior base class and have guassian prior,
// have pointer to prior on this object.
double LangevinDynamics::log_prior(ParameterType &c)
{
    bool single_mode = _opts.single_mode();
    double parameter_sigma = _opts.parameter_proposal_sigma();

    double exponential_constant;
    double variance;
    int    c_dim = parameter_dimension();
    int    M = _opts.cutoff();
    int    idx;
    double log_total = 0;
    // For when we want to infer one fourier mode, the cutoff is M=1
    // which defines 4 modes, but we only want one.
    if((M==1)&&(single_mode))
    {
        // Single mode M=1, final mode is k = (1,1), k^2 = 2.
        variance = pow(4*M_PI*M_PI*2, -1)*parameter_sigma*parameter_sigma;
        log_total -= (0.5/variance)*( pow( std::abs(c(c_dim-1)), 2) );
        return log_total;
    }
    else
    {
        for( size_t i=1; i<M+1; ++i )
        {
            variance = pow(4*M_PI*M_PI*(i*i), -1)*parameter_sigma*parameter_sigma;
            log_total -= (0.5/variance)*pow( std::abs( c(i-1) ), 2); 
        }

        for( int i=-M; i<M+1; ++i )
        for( size_t j=1; j<M+1; ++j )
        { 
            idx = M + (j-1)*(2*M+1) + (i+M);
            variance = pow(4*M_PI*M_PI*(i*i + j*j), -1)*parameter_sigma*parameter_sigma;
            log_total -= (0.5/variance)*pow( std::abs( c(idx) ), 2); 
        }
    }
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

void LangevinDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out )
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
 
    if(m!=M)
    {
        out(k, l, m, 0 ) = xtminus1(0)-gradtminus1(0)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));
        out(k, l, m, 1 ) = xtminus1(1)-gradtminus1(1)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));  
    }
    else
    {
        out(k, l+1, 0, 0 ) = xtminus1(0)-gradtminus1(0)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));
        out(k, l+1, 0, 1 ) = xtminus1(1)-gradtminus1(1)*dt + root_2dt*gsl_ran_gaussian(r, exp( log_sigma ));  
    }
}

void LangevinDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file)
{
    int defining_modes = parameter_dimension(_opts);
    
    // headers
    mcmc_file << "#";
    for( size_t m=0; m<defining_modes; ++m)             
        mcmc_file << m << "\t" << m << "\t";
    mcmc_file << "sigma" << std::endl;
    
    // data
    size_t N = _opts.mcmc_trials();
    size_t burn_in = _opts.burn();
    for( size_t n=burn_in; n<N; ++n )
    {
        for( size_t m=0; m<defining_modes; ++m)
            mcmc_file << std::real(ccc(n,m)) << "\t" << std::imag(ccc(n,m)) << "\t";
        mcmc_file << sss(n) << std::endl;
    }
    mcmc_file.close();
} 
