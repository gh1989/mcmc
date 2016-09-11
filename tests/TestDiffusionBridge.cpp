#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include <string>

#include "BridgeDynamics.h"
#include "LangevinDynamics.h"
#include "Options.h"
#include "Particle.h"

/*
./tests/TestDiffusionBridge -M 30 -Q 500 -o 0.01 -p 0.0002 -K 1 -l 0.1 -P 4
*/

using namespace std;

int main( int argc, char *argv[] )
{
    Options o( argc, argv );
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, o.rng_seed() );
      
    size_t num_particles = o.number_particles();
    
    size_t K = o.parallel_paths();
    size_t L = o.path_length();
    size_t M = o.extra_data_ratio();
   
    std::cout<< "M = " << M << std::endl;
 
    auto ld = LangevinDynamics(o);
    auto  c = ld.default_parameters(o);
    
    if (!o.infer_drift_parameters())
    for(size_t i =0; i<4; ++i)
        c(i) = 0;
    
    std::vector<std::shared_ptr<Particle<BridgeDynamics>>> particles;
    std::vector<std::shared_ptr<Particle<BridgeDynamics>>> resampled;
       
    Particle<BridgeDynamics>::CoarsePathType observed(K, L, 2);
    Particle<BridgeDynamics>::PathType real(K,L,M,2);
   
    double log_sigma = o.log_real_sigma();
    
    Tensor<double, 1> phat(L);
    
    // Randomly dot around points in [-1,1]X[-1,1]
    for(size_t k=0; k<K; ++k)
    {
        real(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        real(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
    
    // Simulate them all forwards.
    for( size_t k=0; k<K; ++k )
    {
        for( size_t l=0; l<L-1; ++l )
        {
        for( size_t m=1; m<M; ++m )
        {
            ld.forward_sim(r, c, log_sigma, observed, k, l, m, real);   
        }
        ld.forward_sim(r, c, log_sigma, observed, k, l+1, 0, real);
        }
    }
    
    // Measure them with noise.
    double observation_sigma = o.observation_noise_sigma();
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
    
    std::cout << "The observations are" << std::endl;

    std::cout << observed;

    Tensor<double, 2> W(L,num_particles);
    Tensor<double, 2> resampled_W(L,num_particles);
    double total = 0;

    for(size_t i=0; i<num_particles; ++i)
    {
        particles.push_back( std::make_shared<Particle<BridgeDynamics>>(o) );
    }
    
    for(size_t i=0; i<num_particles; ++i)
    {
        particles[i]->setup_starts( r, observed );
        W(0, i) = -log(num_particles);
    }

    for(size_t i=0; i<num_particles; ++i) 
        total += exp( W(0, i) ); 
    
    double sigma = o.observation_noise_sigma();
    double a_constant = 0.5 / (sqrt(2*M_PI)*sigma );
    
    for(size_t i=0; i<num_particles; ++i)
        phat(0) += a_constant*exp( particles[i]->unnormal_weight(0, c, log_sigma, observed) );
    phat(0) *= 1.0/num_particles;
    unsigned int resample_particle_index;
    double p[num_particles];
    
    for(size_t i=0; i<num_particles; ++i) 
        p[i] = exp( W(0, i) ) / total; 

    gsl_ran_discrete_t *g;
    g = gsl_ran_discrete_preproc(num_particles, p);
    
    resampled = particles;
    resampled_W = W;
    for( size_t i=0; i<num_particles; ++i )
    {
        resample_particle_index = gsl_ran_discrete(r, g);
        particles[i] = resampled[resample_particle_index];
        W(0,i) = resampled_W(0, resample_particle_index);
        std::cout << "Swapping particle " << i << " for particle " << resample_particle_index << std::endl;   
        std::cout << p[i] <<" against "<< p[resample_particle_index] <<std::endl;
    }

    for(size_t l=1; l<L; ++l)  
    {   
        for(size_t i=0; i<num_particles; ++i)
        {  
            particles[i]->forward_sim(r, c, log_sigma, observed, l); // From [t, t+1)
            W(l,i) = particles[i]->unnormal_weight(l, c, log_sigma, observed);
        } 
        
        // Calculate a part of the marginal.
        total = 0;
        for(size_t i=0; i<num_particles; ++i) 
        {
            std::cout<< "W(" << l<<"," << i << ")" << " = " << W(l,i) << std::endl;
            std::cout<< "exp( W(" << l <<"," << i << "))" << " = " << exp(W(l,i)) << std::endl;
            total += exp( W(l, i) );
        }
        
        resampled = particles;
       
        for( size_t i=0; i<num_particles; ++i )
            p[i] = exp( W(l, i) ) / total;

        gsl_ran_discrete_t *g;
        g = gsl_ran_discrete_preproc(num_particles, p);
        
        resampled_W = W;
        for( size_t i=0; i<num_particles; ++i )
        {
            resample_particle_index = gsl_ran_discrete(r, g);
            particles[i] = resampled[resample_particle_index];
            W(l,i) = resampled_W(l, resample_particle_index);
            std::cout << "Swapping particle " << i << " for particle " << resample_particle_index << std::endl;   
            std::cout << p[i] <<" against "<< p[resample_particle_index] <<std::endl;
        }
    }
    
    std::cout<<"log_marginal "<< phat(L-1) << std::endl;
    
    
    BridgeDynamics::PathType path;

    total = 0;
    for(size_t i=0; i<num_particles; ++i) 
        total += exp( W(0, i) ); 
    for( size_t i=0; i<num_particles; ++i )
        p[i] = exp( W(0, i) ) / total;

    
    for(size_t t=1; t<L; ++t )
    {
        total = 0;
        for(size_t i=0; i<num_particles; ++i) 
            total += exp( W(t, i) ); 
        for( size_t i=0; i<num_particles; ++i )
            p[i] *= exp( W(t, i) ) / total;
    }
    
    /*
     * File output.
     */
     
    std::string ts_file_name = "output/ts_smc_bridge.txt";
    std::ofstream ts_file;
    ts_file.open(ts_file_name);
     
    for(size_t l=0; l<L; ++l)
    {
        // Print for t=l
        for(size_t k=0; k<K; ++k)
            for(size_t i=0; i<num_particles; i++)
            {
                path = particles[i]->path();
                ts_file << path(k,l,0,0) << "\t" << path(k,l,0,1) << "\t";
            }
            
        ts_file << std::endl;
        
        // Fill to t+1
        if(l<L-1)
        for(size_t m=1; m<M; ++m)
        {
            for(size_t k=0; k<K; ++k)
                for(size_t i=0; i<num_particles; i++)
                {
                    path = particles[i]->path();
                    ts_file << path(k,l,m,0) << "\t" << path(k,l,m,1) << "\t";
                }
            ts_file << std::endl;
        }
    }
    ts_file.close();
    
    ts_file_name = "output/ts_smc_bridge_observation.txt";
    ts_file.open(ts_file_name);
    for(size_t l=0; l<L; ++l)
    {
        for(size_t k=0; k<K; ++k)
            ts_file << observed(k,l,0) << "\t" << observed(k,l,1) << "\t";
        
        ts_file<<std::endl;
    }
    ts_file.close();
    
    ts_file_name = "output/ts_smc_bridge_probs.txt";
    ts_file.open(ts_file_name);
    for(size_t i=0;i<num_particles;++i)
    {
        ts_file << p[i] << "\t";
    }
    ts_file << std::endl;
    ts_file.close();
    
    particles.clear();
    resampled.clear();

}
