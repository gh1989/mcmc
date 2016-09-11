#include "CurvedSurfaceDynamics.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double CurvedSurfaceDynamics::log_p( CoarsePathType &y,  PathType &x, ParameterType &c, double log_sigma, int k, int l)
{
    
    double obs_sigma =  _opts.observation_noise_sigma();
    double exponential_constant_obs  = 0.5/(obs_sigma*obs_sigma);
    double log_total = 0;

    log_total -= exponential_constant_obs * pow( x(k, l, 0, 0 ) - y(k, l, 0), 2);
    log_total -= exponential_constant_obs * pow( x(k, l, 0, 1 ) - y(k, l, 1), 2);
    
    return log_total;
}

ComplexType CurvedSurfaceDynamics::sample_transition_density(gsl_rng *r, ComplexType c )
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    double random_real_part = gsl_ran_gaussian(r, proposal_sigma);
    double random_imag_part = gsl_ran_gaussian(r, proposal_sigma);
    return ComplexType(random_real_part, random_imag_part) + c;
}

CurvedSurfaceDynamics::ParameterType CurvedSurfaceDynamics::sample_transition_density(gsl_rng *r, ParameterType& C)
{
    double proposal_sigma = _opts.parameter_proposal_sigma();
    bool single_mode = _opts.single_mode();
    
    size_t c_dim = parameter_dimension(_opts);
    ParameterType C_star(c_dim);
    
    for(size_t i=0; i<c_dim; ++i) 
        C_star(i) = 0.0;
    
    int M = _opts.cutoff();
    int idx;
    double s = 2;
    
    double eigenvalue_power;
    if ((M==1)&&(single_mode))
    {
        eigenvalue_power = 2*M_PI*sqrt(2);
        C_star(c_dim-1) = ComplexType(  gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ), 
                                        gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ) ) + C(c_dim-1); 
    }
    else
    {
        for( size_t i=1; i<M+1; ++i )
        {
            eigenvalue_power = 2*M_PI*(i);
            C_star(i-1) = ComplexType( gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ), 
                                       gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ) ) + C(i-1); 
        }

        for( int i=-M; i<M+1; ++i )
            for( size_t j=1; j<M+1; ++j )
            { 
                idx = M + (j-1)*(2*M+1) + (i+M);
                eigenvalue_power = 2*M_PI*sqrt(i*i + j*j);
                C_star(idx) = ComplexType( gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ), 
                                           gsl_ran_gaussian(r, proposal_sigma*pow(eigenvalue_power, -s) ) ) + C(idx); 
            }
    }
    
    if (gsl_rng_uniform(r)<0.5)
        C_star=-C_star;  
    
    return C_star;
}

double CurvedSurfaceDynamics::log_transition(ParameterType &c_star, ParameterType &c)
{
    return 0;
}

double CurvedSurfaceDynamics::log_prior(ParameterType &c)
{
    bool single_mode = _opts.single_mode();
    double parameter_sigma = _opts.parameter_proposal_sigma();

    double exponential_constant;
    double variance;
    int    c_dim = parameter_dimension();
    int    M = _opts.cutoff();
    int    idx;
    double log_total = 0;
    
    double s = 2;
    
    if((M==1)&&(single_mode))
    {
        // Single mode M=1, final mode is k = (1,1), k^2 = 2.
        // c ~ N(0, (4 pi^2 k^2)^{ -s }.
        
        variance = pow(4*M_PI*M_PI*2, -s)*parameter_sigma*parameter_sigma;
        log_total -= (0.5/variance)*( pow( std::abs(c(c_dim-1)), 2) );
        return log_total;
    }
    else
    {
        for( size_t i=1; i<M+1; ++i )
        {
            variance = pow(4*M_PI*M_PI*(i*i), -s)*parameter_sigma*parameter_sigma;
            log_total -= (0.5/variance)*pow( std::abs( c(i-1) ), 2); 
        }

        for( int i=-M; i<M+1; ++i )
        for( size_t j=1; j<M+1; ++j )
        { 
            idx = M + (j-1)*(2*M+1) + (i+M);
            variance = pow(4*M_PI*M_PI*(i*i + j*j), -s)*parameter_sigma*parameter_sigma;
            log_total -= (0.5/variance)*pow( std::abs( c(idx) ), 2); 
        }
    }
    return log_total;
}

double CurvedSurfaceDynamics::log_path_likelihood(  PathType &x, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m )
{
    double log_total = 0;
    double dt = _opts.trajectory_path_delta();
    double obs_sigma = _opts.observation_noise_sigma();

    Vector2d xt;
    Vector2d xtminus1;
    Vector2d drift_value;

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
        double exponential_constant_path = 0.25 / ( pow(diffusion(xtminus1,c),2)*dt );
        drift_value = drift( xtminus1, c);
        log_total -= exponential_constant_path * pow( xt(0)-xtminus1(0)-drift_value(0)*dt, 2);
        log_total -= exponential_constant_path * pow( xt(1)-xtminus1(1)-drift_value(1)*dt, 2);
    }
    
    return log_total;
}

Vector2d CurvedSurfaceDynamics::drift( Vector2d &x, ParameterType &c )
{
    FourierSeries h(_cutoff);
    h.set_modes( c );
    
    // TODO: optimise this function call.
    auto dh = h.grad( x );
   
    double h_x = dh(0);
    double h_y = dh(1);
    
    // TODO: optimise this function call.
    auto hessian_h = h.hessian( x );
    
    double h_xx = hessian_h(0,0);
    double h_yy = hessian_h(1,1);
    double h_xy = hessian_h(0,1);  
    
    Vector2d grad_det_g;
    
    grad_det_g(0) = -h_xy * pow(h_x,3) + 2*h_y*h_xx*pow(h_x,2) + 2*pow(h_y,2)*h_xx*h_x+2*h_y*h_xy;
    grad_det_g(1) = -h_yy * pow(h_x,3) + 2*h_y*h_xy*pow(h_x,2) + 2*pow(h_y,2)*h_xy*h_x+2*h_y*h_yy;
    
    // I copy and pasted this from somewhere else in the code (*)
    double det_g = 1+pow(h_y,2)+pow(h_x,2) + pow(h_y,2)*pow(h_x,2)-pow(h_x,3)*h_y;
   
    Matrix2d A(2,2);
    A << 1 + pow( h_y, 2), -h_x*h_y, -h_y*h_x, 1 + pow( h_x, 2);
    
    Vector2d divA;
    divA << 2*h_y*h_xy - h_xy*h_y-h_y*h_xx, 2*h_x*h_xy - h_xy*h_y-h_x*h_yy;
    
    double sigma = exp( _opts.log_real_sigma() );
    
    return sigma*( pow( det_g, -0.5 ) * divA - 0.5 * pow( det_g, -1.5 )*(A*grad_det_g) );
}

double CurvedSurfaceDynamics::diffusion( Vector2d &x, ParameterType &c )
{
    FourierSeries h(_cutoff);
    h.set_modes(c);
    
    // TODO: optimise this function call.
    auto dh = h.grad( x );
    double h_x = dh(0);
    double h_y = dh(1);
    
    // TODO: optimise this function call.
    auto hessian_h = h.hessian( x );
    
    double h_xx = hessian_h(0,0);
    double h_yy = hessian_h(1,1);
    double h_xy = hessian_h(0,1);  
    
    // I copy and pasted this from somewhere else in the code (*)
    double det_g = 1+pow(h_y,2)+pow(h_x,2) + pow(h_y,2)*pow(h_x,2)-pow(h_x,3)*h_y;
    
    double sigma = exp( _opts.log_real_sigma() );
    
    return sigma*det_g;
}

void CurvedSurfaceDynamics::forward_sim( gsl_rng *r, ParameterType &c, double log_sigma, CoarsePathType &y, int k, int l, int m, PathType &out )
{
    int    M = _opts.extra_data_ratio();
    double dt = _opts.trajectory_path_delta();
    double root_2dt = sqrt(2.0*dt);

    FourierSeries V(_cutoff);
    V.set_modes(c);
            
    Vector2d xtminus1;
    Vector2d drift_value = drift( xtminus1, c );
            
    // If (m==0) then take the value from the last observation interval.
    if( m > 0)
        xtminus1 << out(k, l, m-1, 0 ), out(k, l, m-1, 1 );
    else
        xtminus1 << out(k, l-1, M-1, 0), out(k, l-1, M-1, 1);
 
    double diffusion_coefficient = diffusion( xtminus1, c );
 
    if (m==M)
    {
        m = 0;
        l += 1;
    }
    
    // Forward simulation.
    out(k, l, m, 0 ) = xtminus1(0)+drift_value(0)*dt + root_2dt*gsl_ran_gaussian(r, diffusion_coefficient);
    out(k, l, m, 1 ) = xtminus1(1)+drift_value(1)*dt + root_2dt*gsl_ran_gaussian(r, diffusion_coefficient);  
    
    // Impose the boundary conditions on [0,1]X[0,1].
    out(k, l, m, 0 ) -= int(  out(k, l, m, 0 ) );
    out(k, l, m, 1 ) -= int(  out(k, l, m, 1 ) ); 
    
}

void CurvedSurfaceDynamics::output_file_timeseries(ParameterChainType &ccc, SigmaChainType &sss, std::ofstream &mcmc_file)
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
