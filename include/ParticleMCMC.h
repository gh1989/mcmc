#ifndef PARTICLE_MCMC_H
#define PARTICLE_MCMC_H

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
    
    double smc(gsl_rng *r);    
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
    double log_marginal_likelihood_c_star;
    double log_marginal_likelihood_c;
    
    Tensor<double, 2> W;
  
    std::vector<std::shared_ptr<Particle<Dynamics_>>> particles;
    std::vector<std::shared_ptr<Particle<Dynamics_>>> resampled;
};
      
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::accept() 
{  
    c = c_star;
    log_sigma = log_sigma_star;
    log_marginal_likelihood_c = log_marginal_likelihood_c_star;
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::finish() 
{
    _dynamics.output_file_timeseries( parameter_chain, log_sigma_chain );
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::generate_observations(gsl_rng *r)
{   
    double observation_sigma = _opts.observation_noise_sigma();
    double random_noise_x;
    double random_noise_y;

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
    
    log_marginal_likelihood_c = smc(r);
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::generate_random_starts(gsl_rng *r, PathType &out )
{
    for( size_t k=0; k<K; ++k )
    {
        // Initialisation: random point inside [-1,1]X[-1,1].
        out(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        out(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
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
    if( infer_drift_parameters )
    {
        for( size_t i=0; i<c_dim; ++i)
            c_star(i) = _dynamics.sample_transition_density(r, c(i));
    }
    else
        c_star = c;
    
    double log_sigma_proposal_standard_deviation = _opts.parameter_proposal_diffusion_sigma();
    
    if( infer_diffusion_parameters )
        log_sigma_star = gsl_ran_gaussian(r, log_sigma_proposal_standard_deviation );
    else
        log_sigma_star = log_sigma;

    log_marginal_likelihood_c_star = smc( r );
    std::cout << "log_marginal_likelihood_c_star: " << log_marginal_likelihood_c_star << std::endl;
    std::cout << "log_marginal_likelihood_c: " << log_marginal_likelihood_c << std::endl;
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::setup_observed_starts(  gsl_rng *r, 
                                                      CoarsePathType &y, 
                                                      PathType &out )
{
    double sigma = _opts.observation_noise_sigma();
    //std::cout << "ParticleMCMC thinks that observation_noise_sigma is: " << sigma << std::endl;
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
          _dynamics.forward_sim(r, c_, log_sigma_, k, l, m, out);   
        }
        _dynamics.forward_sim(r, c_, log_sigma_, k, l+1, 0, out);
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
    double log_total = 0;    
    // Sigma
    if( infer_diffusion_parameters )
    {
        log_total += _dynamics.log_transition_sigma(log_sigma, log_sigma_star);
        log_total -= _dynamics.log_transition_sigma(log_sigma_star, log_sigma);
        log_total += _dynamics.log_prior_sigma(log_sigma_star);
        log_total -= _dynamics.log_prior_sigma(log_sigma);
    }
    
    // symmetric
    // log_total += _dynamics.log_transition(c, c_star);       
    // log_total -= _dynamics.log_transition(c_star, c);
    log_total += _dynamics.log_prior(c_star);
    log_total -= _dynamics.log_prior(c); 
        
    log_total += log_marginal_likelihood_c_star;
    log_total -= log_marginal_likelihood_c;
    
    return log_total;
    
}

template<class Dynamics_>
double ParticleMCMC<Dynamics_>::smc(gsl_rng *r)
{
    double total;
    Tensor<double, 1> phat(L);
    
    for(size_t i=0; i<num_particles; ++i)     
        particles[i]->setup_starts( r, observed );
    
    for(size_t i=0; i<num_particles; ++i)  
        W(0, i) = -log(num_particles);
        
    total = 0;
    for(size_t i=0; i<num_particles; ++i) 
        total += exp( W(0, i) ); 
       
    phat(0) = (1.0/num_particles)*total;
    resample_with_replacement(r, 0);   
    
    std::cout << "c_star" << c_star << std::endl;
    std::cout << "log_sigma_star" << log_sigma_star << std::endl;

    for(size_t t=1; t<L; ++t)
    {
        // Generate particle samples.
        for(size_t i=0; i<num_particles; ++i)     
        {
            particles[i]->forward_sim(r, c_star, log_sigma_star, t );
            W(t, i) = particles[i]->unnormal_weight(t, c_star, log_sigma_star, observed);
        }
        
        // Calculate a part of the marginal.
        total = 0;
        for(size_t i=0; i<num_particles; ++i) 
            total += exp( W(t, i) );
        
        phat(t) = phat(t-1) * (1.0/num_particles)*total;
        
        // Resample the particles based on weights.
        resample_with_replacement(r, t);
    }

    // Return the marginal estimate.
    return log( phat(L-1) );
}

template<class Dynamics_>
void ParticleMCMC<Dynamics_>::resample_with_replacement(gsl_rng *r, size_t t)
{    
    resampled = std::move( particles ); // particles of size zero now.
    unsigned int resample_particle_index;
    double p[num_particles];
    for( size_t i=0; i<num_particles; ++i )
        p[i] = exp( W(t, i) );

    gsl_ran_discrete_t *g;
    g = gsl_ran_discrete_preproc(num_particles, p);
   
    for( size_t i=0; i<num_particles; ++i )
    {
        resample_particle_index = gsl_ran_discrete(r, g);
        particles.push_back( resampled[resample_particle_index] );
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
