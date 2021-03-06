#ifndef PARTICLE_MCMC_H
#define PARTICLE_MCMC_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

#include "Particle.h"

using namespace MCMC;

namespace MCMC
{

template< class Dynamics_>
class ParticleMCMC
{
public:
    using ParameterType         = typename Dynamics_::ParameterType;
    using ParameterChainType    = typename Dynamics_::ParameterChainType;
    using PathType              = typename Dynamics_::PathType;
    using CoarsePathType        = typename Dynamics_::CoarsePathType;
    using SigmaChainType        = typename Dynamics_::SigmaChainType;    

    ParticleMCMC(Options& o) : _dynamics(o)
    {
        setup_from_options(o);
    } 
    
    // Do not instantiate without options.
    ParticleMCMC()=delete;
    
    // No need for this 
    ParticleMCMC(ParticleMCMC& other) : _dynamics(other.dynamics()) 
    {
        std::cout<<"ParticleMCMC copy constructor called in " << this <<std::endl;
        setup_from_options( other.opts() );
    }
    
    ParticleMCMC& operator=(ParticleMCMC& other)    
    {
        std::cout<<"ParticleMCMC copy assignment operator= called in " << this <<std::endl;
        _dynamics = other.dynamics();
        setup_from_options( other.opts() );
    }
        
    ParticleMCMC(ParticleMCMC&& other) : _dynamics(other.dynamics())
    {
        std::cout<<"ParticleMCMC move constructor called in " << this <<std::endl;
        setup_from_options( other.opts() );
    }
    
    ParticleMCMC& operator=(ParticleMCMC&& other)
    {
        std::cout<<"ParticleMCMC move assignment operator= called in " << this <<std::endl;
        _dynamics = other.dynamics();
        setup_from_options(other.opts());
        return *this;
    };
    
    ~ParticleMCMC()
    {
        std::cout<<"~ParticleMCMC() called for: "<< this << std::endl;
    };
    
    void accept();
    void finish();
    void generate_observations( gsl_rng *r);
    void generate_true_trajectory(gsl_rng *r);
    void propose(gsl_rng *r); 
    void store_chain(int n);
    void resample_with_replacement(gsl_rng *r, size_t t);
    void setup_from_options(Options &o);
    void log_particle_end_points( std::ofstream &endpoints_file);
    
    ParameterType& current_drift(){ return c; }
    double current_log_sigma() { return log_sigma; }
    
    long double smc(gsl_rng *r);    
    double log_acceptance_probability();
    
    Options& opts(){ return _opts; }
    Dynamics_& dynamics(){ return _dynamics;}
    
private:
    void generate_random_starts( gsl_rng *r, PathType &out );
    void setup_observed_starts( gsl_rng *r, CoarsePathType &y, PathType &out );
    void trajectory( gsl_rng *r, ParameterType &c, double sigma_, PathType &out );
    
    size_t num_particles;

    Options                 _opts;   
    Dynamics_               _dynamics;
    
    size_t K;
    size_t L;
    size_t M;
    size_t N;

    CoarsePathType          observed;
    PathType                real;
    PathType                x;
    PathType                x_star;
    
    size_t                  c_dim;
    ParameterType           c;
    ParameterType           c_star;
    ParameterType           real_c;
    ParameterChainType      parameter_chain;
    
    double                  log_real_sigma;
    double                  log_sigma_star;
    double                  log_sigma;  
    SigmaChainType          log_sigma_chain; 
    
    bool    infer_drift_parameters;
    bool    infer_diffusion_parameters;
  
    // SMC/Particle vars
    long double marginal_likelihood_c_star;
    long double marginal_likelihood_c;
    
    Tensor<double, 2> W;
     
    std::ofstream mcmc_file; 
    
    std::vector<std::shared_ptr<Particle<Dynamics_>>> particles;
    std::vector<std::shared_ptr<Particle<Dynamics_>>> resampled;
};
      
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::accept() 
{  
    c = c_star;
    log_sigma = log_sigma_star;
    marginal_likelihood_c = marginal_likelihood_c_star;
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::finish() 
{
    
    unsigned long int milliseconds_since_epoch = 
    std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
        
    std::string filename = "output/";
    filename += _opts.output_subfolder();
    filename += "/pMCMCTimeSeries_";
    filename += Dynamics_::dynamics_string();
    filename += std::to_string( milliseconds_since_epoch );
    filename += ".txt";
    
    mcmc_file.open(filename);
    _opts.print_header( mcmc_file );
    _dynamics.output_file_timeseries( parameter_chain, log_sigma_chain, mcmc_file );
}


template<class Dynamics_>
void ParticleMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
    double observation_sigma = _opts.observation_noise_sigma();
    
    double random_noise_x;
    double random_noise_y;

    // Generate observations from the real trajectories.
    for( size_t k=0; k<K; ++k )
    {
        for( size_t l=0; l<L; ++l )
        {   
            random_noise_y = gsl_ran_gaussian(r, observation_sigma );
            random_noise_x = gsl_ran_gaussian(r, observation_sigma );          
            
            observed(k, l, 0 ) = real(k, l, 0, 0 ) + random_noise_x;           
            observed(k, l, 1 ) = real(k, l, 0, 1 ) + random_noise_y;
        }   
    }
    
    // Run SMC and store marginal_likelihood for time 0.
    marginal_likelihood_c = smc(r);
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [0,1]X[0,1].
        out(k, 0, 0, 0 ) = gsl_rng_uniform(r);
        out(k, 0, 0, 1 ) = gsl_rng_uniform(r);
    }
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::log_particle_end_points( std::ofstream &endpoints_file )
{
    // A utility function for logging the ends of trajectories.
    PathType the_path;
    for( int i=0; i<num_particles; i++)
    {
        the_path = particles[i]->path(); 
        for(int k=0; k<K; ++k)
            endpoints_file << the_path(k,L-1,0,0) << "\t" << the_path(k,L-1,0,1) << "\t";
    }
    endpoints_file << std::endl;
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::generate_true_trajectory(gsl_rng *r)
{ 
    generate_random_starts( r, real );
    trajectory(r, real_c, log_real_sigma, real );
}


template<class Dynamics_>
void ParticleMCMC<Dynamics_>::propose( gsl_rng *r )
{
    // With what standard deviation to propose log(sigma_d)?
    double log_sigma_proposal_standard_deviation = _opts.parameter_proposal_diffusion_sigma();
    
    // If inferring drift parameters (which actually means inferring modes, in the curved surface implementation
    // these modes are included in the drift term, too.
    c_star = infer_drift_parameters ? _dynamics.sample_transition_density(r, c) : c;    
    
    // If inferring diffusion constant...
    log_sigma_star = infer_diffusion_parameters ? (gsl_ran_gaussian(r, log_sigma_proposal_standard_deviation ) + log_sigma) : log_sigma;
    
    // Setup the proposal in line with observations for parallel chains.
    setup_observed_starts( r, observed, x_star );
    
    // Simulate trajectories.
    trajectory( r, c_star, log_sigma_star, x_star );

    // Run SMC store new marginal likelihood.
    marginal_likelihood_c_star = smc( r );    
    
    std::cout << "[Debug] marginal_likelihood_c_star: " << marginal_likelihood_c_star   << std::endl;
    std::cout << "[Debug] marginal_likelihood_c: "      << marginal_likelihood_c        << std::endl;       
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::setup_observed_starts(  gsl_rng *r, 
                                                      CoarsePathType &y, 
                                                      PathType &out )
{
    double sigma = _opts.observation_noise_sigma();
    for( size_t k=0; k<K; ++k )
    {
        out(k, 0, 0, 0 ) = y(k, 0, 0);
        out(k, 0, 0, 1 ) = y(k, 0, 1);
    }
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::trajectory(   gsl_rng *r, 
                                            ParameterType &c_,
                                            double log_sigma_,
                                            PathType &out )
{
    for( size_t k=0; k<K; ++k )
      for( size_t l=0; l<L-1; ++l )
      {
        for( size_t m=1; m<M; ++m )
        {
          _dynamics.forward_sim(r, c_, log_sigma_, observed, k, l, m, out);   
        }
        _dynamics.forward_sim(r, c_, log_sigma_, observed, k, l+1, 0, out);
      }
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::store_chain(int n)
{
    for( size_t i=0; i<c_dim; ++i)
        parameter_chain(n, i) = c(i);
    log_sigma_chain(n) = log_sigma;
}

template<class Dynamics_>
double ParticleMCMC<Dynamics_>::log_acceptance_probability()
{
    long double log_total = 0;    

    if( infer_diffusion_parameters )
    {
        log_total += _dynamics.log_prior_sigma(log_sigma_star);
        log_total -= _dynamics.log_prior_sigma(log_sigma);
    }
    
    if( infer_drift_parameters )
    {
        log_total += _dynamics.log_prior(c_star);
        log_total -= _dynamics.log_prior(c); 
    }
   
    log_total += log( marginal_likelihood_c_star / marginal_likelihood_c );
    
    return log_total;
    
}

template<class Dynamics_>
long double ParticleMCMC<Dynamics_>::smc(gsl_rng *r)
{
    double sigma = _opts.observation_noise_sigma();
    double a_constant = 0.5 / ( pow(2*M_PI*sigma*sigma, 1.5) ); // WHY was this not 3/2 ???
    double current_mean;
    double total;
    
    // The marginal likelihood estimate evolves in time, updated after each SMC loop below.
    Tensor<long double, 1> marginal_likelihood_estimate(L);
    
    // Align all particles at the start of the observed.
    for(size_t i=0; i<num_particles; ++i)     
        particles[i]->setup_starts( r, observed );
    
    // Calculate the unnormalised weights for each partical, store in W(0,i).
    for(size_t i=0; i<num_particles; ++i)  
        W(0, i) = particles[i]->unnormal_weight(0, c_star, log_sigma_star, observed);
        
    // Calculate the mean of the unnormalised weights at time 0.
    total = 0;
    for(size_t i=0; i<num_particles; ++i) 
        total += a_constant*exp( W(0, i) );   
    current_mean = total / num_particles;
    
    // Marginal likelihood estimate at time 0.
    marginal_likelihood_estimate(0) = current_mean;
    
    // Resample.
    resample_with_replacement(r, 0);   

    for(size_t t=1; t<L; ++t)
    {
        // Generate particle samples.
        for(size_t i=0; i<num_particles; ++i)     
        {
            // Forward simulate particle i, from t-1 to t incrementally in M steps.
            particles[i]->forward_sim(r, c_star, log_sigma_star, observed, t );
            
            // The unnormalised log weight at time t for particle i. p( y | x, c ).
            W(t, i) = particles[i]->unnormal_weight(t, c_star, log_sigma_star, observed);
        }
        
        // Calculate a part of the marginal.
        total = 0;
        for(size_t i=0; i<num_particles; ++i) 
            total += a_constant*exp( W(t, i) );   
        current_mean = total / num_particles;
        
        // \hat{p}(y(t)|x) =\hat{p}(y(t-1)|x) * mean( unnormalised weights ).
        marginal_likelihood_estimate(t) = marginal_likelihood_estimate(t-1)*current_mean;   
        
        resample_with_replacement(r, t);
    }

    // Return the log marginal estimate.
    return marginal_likelihood_estimate(L-1);
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::resample_with_replacement(gsl_rng *r, size_t t)
{    
    // Move particles to resampled.
    resampled = std::move( particles );
    
    unsigned int resample_particle_index;
    double p[num_particles];
    
    // Set of unnormalised weights.
    for( size_t i=0; i<num_particles; ++i )
        p[i] = exp( W(t, i) );

    gsl_ran_discrete_t *g;
    g = gsl_ran_discrete_preproc(num_particles, p);
   
    // Create a copy of the weights
    auto resampled_W = W;
    
    for( size_t i=0; i<num_particles; ++i )
    {
        // Random index to resample a particle
        resample_particle_index = gsl_ran_discrete(r, g);
        particles.push_back( resampled[resample_particle_index] );
        
        // Push back puts it on the end of the vector...
        W(t, num_particles - i - 1) = resampled_W(t, resample_particle_index);
    }
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::setup_from_options(Options &o)
    {
        _opts = o;

        // The number of parallel trajectories per iteration.
        K = o.parallel_paths();
        
        // The length of each trajectory path in timesteps.
        L = o.path_length();
        
        // The extra data in between observed data ratio.
        M = o.extra_data_ratio();
        
        // The number of iterations for the MCMC experiment.
        N = o.mcmc_trials();

        // The observations K*L*2 Tensor of doubles
        observed    = CoarsePathType( K, L, 2 );
        
        // The actual path
        real        = PathType( K, L, M, 2 );
        
        // The proposed and current path for the MCMC simulation
        x           = PathType( K, L, M, 2 );
        x_star      = PathType( K, L, M, 2 );    
        
        // We may want to infer these separately
        infer_drift_parameters = o.infer_drift_parameters();
        infer_diffusion_parameters = o.infer_diffusion_parameters();
        
        // The trajectories: dX_t = b(X_t) dt + sqrt(2*D(X_t)) dW_t
        // b(X_t) : related parameters
        
        // The dimension of the b(X_t) parameters
        c_dim               = _dynamics.parameter_dimension(o);

        // The current parameters for b(X_t)
        c                   = ParameterType(c_dim);
        // The proposed parameters for b(X_t)
        c_star              = ParameterType(c_dim);
        
        // The genuine parameters that need to be inferred
        real_c              = _dynamics.default_parameters(o);
        
        if( infer_drift_parameters )
            for(size_t i=0; i<c_dim; ++i)  c(i) = 0;
        else
            c = real_c;
                
        // The chain of b(X_t) parameters
        parameter_chain     = ParameterChainType(N, c_dim);
            
        // D(X_t) : in this case just the diffusion coefficient, distributed log normally so we
        // propose log_sigma with Gaussian prior
        log_real_sigma              = o.log_real_sigma();
        
        if( infer_diffusion_parameters )
        {
            std::cout<<"Inferring diffusion parameters."<<std::endl;
            log_sigma               = o.log_start_sigma();
        }
        else
            log_sigma = log_real_sigma;
        
        log_sigma_chain         = SigmaChainType(N);    
        
        // SMC/Particle vars
        num_particles = o.number_particles();
        W = Tensor<double, 2>(L, num_particles);
        for(size_t i=0; i<num_particles; ++i)
        {
            particles.push_back( std::make_shared<Particle<Dynamics_>>(o) );
        }
        

      
    }

#endif
