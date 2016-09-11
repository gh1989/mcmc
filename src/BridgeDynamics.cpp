#include "BridgeDynamics.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


double BridgeDynamics::log_p( CoarsePathType &y,  PathType &x, ParameterType &c, double log_sigma, int k, int l)
{
    
    double obs_sigma =  _opts.observation_noise_sigma();

    double exponential_constant_obs  = 0.5/(obs_sigma*obs_sigma);
    double log_total = 0;

    log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
    log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    
    return log_total;
}


ComplexType BridgeDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    double random_real_part = gsl_ran_gaussian(r, proposal_sigma);
    double random_imag_part = gsl_ran_gaussian(r, proposal_sigma);
    return ComplexType(random_real_part, random_imag_part) + c;
}

BridgeDynamics::ParameterType BridgeDynamics::sample_transition_density(gsl_rng *r, ParameterType& C)
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

double BridgeDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    return 0;
}

double BridgeDynamics::log_prior(ParameterType &c)
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

// Poorly named: this is just the probability of p( x_{l*M*dt+m*dt} | x_{l*M*dt+m*dt - dt}, c, log_sigma, y )
// unless m==0 then its just the conditional prob p( x_{l*M*dt} | y_{l*M*dt}, c )
double BridgeDynamics::log_path_likelihood(  PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m )
{
    double log_total = 0;
    double dt = _opts.trajectory_path_delta();
    size_t M = _opts.extra_data_ratio();
    double obs_sigma = _opts.observation_noise_sigma();
    FourierSeries V(_cutoff);
    V.set_modes( c );

    Vector2d xt;
    Vector2d xtminus1;
    Vector2d gradvtminus1;

    double exponential_constant_path = 0.25 / ( exp(2*log_sigma)*dt );
    double exponential_constant_obs = 0.5 / ( obs_sigma*obs_sigma );
    
    if( m == 0 )
    {        
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
        log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    }
    else
    {
        if( m==M )
            xt << x(k, l+1, 0, 0 ), x(k, l+1, 0, 1 );
        else
            xt << x(k, l, m, 0 ), x(k, l, m, 1 );
        
        xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
        
        gradvtminus1 = V.grad( xtminus1 );
        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(0)+gradvtminus1(0)*dt, 2.0);
        log_total -= exponential_constant_path * pow( xt(1)-xtminus1(1)+gradvtminus1(1)*dt, 2.0);
    }
    
    if( log_total ==0 )
    {
        std::cout<< "WARNING. log_total is 0." << std::endl;
    }
  
    return log_total;
}

// Poor name, this is the probability p^( x^k_{l*M*dt+m*dt} | x^k_{l*M*dt+m*dt - dt}, c, log_sigma, y )
double BridgeDynamics::log_bridge_likelihood( PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m)
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);
    
    FourierSeries V(_cutoff);
    V.set_modes(c);
    
    double bridge_variance;
    // M*dt - m*dt??
    double delta_m = 1.0 - m*dt;
    //double delta_m =  M*dt - m*dt;
    double o = _opts.observation_noise_sigma();
    double s = exp( log_sigma );
    double sigma2tilde = pow( s, 2.0 ) / ( pow( s, 2.0 )*delta_m + pow( o, 2.0 )  );
    
    Vector2d xtminus1;
    Vector2d drift_vector;
    Vector2d end_point;
    Vector2d gradtminus1;
    Vector2d current;
    
    if( m == M )
        current << x(k,l+1,0,0), x(k,l+1,0,0);
    else
        current << x(k,l,m,0), x(k,l,m,0);    
    
    xtminus1 << x(k, l, m-1, 0 ), x(k, l, m-1, 1 );
    end_point << y(k,l+1,0), y(k,l+1,1);
    
    std::cout<< "end_point" << end_point << std::endl;
    std::cout<< "current" << current << std::endl;
    std::cout<< "xtminus1" << xtminus1 << std::endl;
    std::cout<< "sigma2tilde" << sigma2tilde << std::endl;
    
    gradtminus1 = V.grad( xtminus1 );
    std::cout<< "gradtminus1" << gradtminus1 << std::endl;
    std::cout<< "delta_m" << delta_m << std::endl;
    
    drift_vector = -gradtminus1 + sigma2tilde*( end_point - xtminus1 + gradtminus1 * delta_m );
    bridge_variance = pow( s, 2.0 ) * (1 - sigma2tilde*dt )*dt;
 
    double log_total = 0;
    double constant = 0.25 / bridge_variance;
    
    std::cout<< "bridge_variance" << bridge_variance << std::endl;
    
    log_total -= constant*pow( current(0) - xtminus1(0) - drift_vector(0)*dt, 2.0 );
    log_total -= constant*pow( current(1) - xtminus1(1) - drift_vector(1)*dt, 2.0 );
      
    return log_total;
    
}


void BridgeDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out )
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);
    
    FourierSeries V(_cutoff);
    V.set_modes(c);
    
    double bridge_sigma;
    // M*dt - m*dt??
    double delta_m = 1.0 - m*dt;
    //double delta_m =  M*dt - m*dt;
    double o = _opts.observation_noise_sigma();
    double s = exp( log_sigma );
    double sigma2tilde = pow( s, 2.0 ) / ( pow( s, 2.0 )*delta_m + pow( o, 2.0 ) );
    
    Vector2d xtminus1;
    Vector2d drift_vector;
    Vector2d end_point;
    Vector2d gradtminus1;

    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);
    
    end_point << y(k,l+1,0), y(k,l+1,1);
    
    gradtminus1 = V.grad( xtminus1 );
    drift_vector = -gradtminus1 + sigma2tilde*( end_point - xtminus1 + gradtminus1 * delta_m );
    bridge_sigma = root_2dt*s*sqrt(1.0 - sigma2tilde*dt);
    
    if( sqrt(drift_vector.transpose()*drift_vector) > pow(10,10) )
    {        
        std::cout<< "WARNING... large drift?" << std::endl;
        std::cout<< "m " << m << std::endl;
        std::cout<< "l " << l << std::endl;
        std::cout<< "drift_vector" << drift_vector << std::endl;
        std::cout<< "delta_m" << delta_m << std::endl;
        std::cout<< "xtminus1" << xtminus1 << std::endl;
        std::cout<< "gradtminus1" << gradtminus1 << std::endl;
        std::cout<< "end_point"     << end_point << std::endl;
        std::cout<< "sigma2tilde" << sigma2tilde << std::endl;
        std::cout<< "bridge_sigma" << bridge_sigma << std::endl;
        
    }
    
    if(m!=M)
    {
        out(k, l, m, 0 ) = xtminus1(0) + drift_vector(0)*dt + gsl_ran_gaussian(r, bridge_sigma);
        out(k, l, m, 1 ) = xtminus1(1) + drift_vector(1)*dt + gsl_ran_gaussian(r, bridge_sigma);  
    }
    else
    {
        out(k, l+1, 0, 0 ) = xtminus1(0) + drift_vector(0)*dt + gsl_ran_gaussian(r, bridge_sigma);
        out(k, l+1, 0, 1 ) = xtminus1(1) + drift_vector(1)*dt + gsl_ran_gaussian(r, bridge_sigma);      
    }
}

void BridgeDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file)
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
