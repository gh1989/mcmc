#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include <string>

#include "LangevinDynamics.h"
#include "Options.h"
#include "Particle.h"

using namespace std;

/*
The purpose of these tests is to work out what is breaking
*/

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
    
    auto ld = LangevinDynamics(o);
    auto  c = ld.default_parameters(o);
    
    std::vector<std::shared_ptr<Particle<LangevinDynamics>>> particles;
    std::vector<std::shared_ptr<Particle<LangevinDynamics>>> resampled;
       
    Particle<LangevinDynamics>::CoarsePathType observed(K, L, 2);
    Particle<LangevinDynamics>::PathType real(K,L,M,2);
    
    double log_sigma = -6;
    size_t t = 1;
       
    // Randomly dot around points in [-1,1]X[-1,1]
    for(size_t k=0; k<K; ++k)
    {
        real(k, 0, 0, 0 ) = 2.0*gsl_rng_uniform(r)-1;
        real(k, 0, 0, 1 ) = 2.0*gsl_rng_uniform(r)-1;
    }
    
    // Simulate them all forwards.
    for( size_t k=0; k<K; ++k )
        for( size_t l=0; l<L-1; ++l )
        {
        for( size_t m=1; m<M; ++m )
        {
            ld.forward_sim(r, c, log_sigma, k, l, m, real);   
        }
        ld.forward_sim(r, c, log_sigma, k, l+1, 0, real);
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
    
    Tensor<double, 2> W(L,num_particles);
    double total = 0;

    for(size_t i=0; i<num_particles; ++i)
    {
        particles.push_back( std::make_shared<Particle<LangevinDynamics>>(o) );
    }
    
    for(size_t i=0; i<num_particles; ++i)
    {
        particles[i]->setup_starts( r, observed );
        //W(0, i) = -log(num_particles);
    }

    for(size_t i=0; i<num_particles; ++i) 
        total += exp( W(0, i) ); 
    
    unsigned int resample_particle_index;
    double p[num_particles];
    
    for(size_t t=1; t<L; ++t)  
    {   
        for(size_t i=0; i<num_particles; ++i)
        {  
            particles[i]->forward_sim(r, c, log_sigma, t); // It was here... t changed to t-1
            W(t,i) = particles[i]->unnormal_weight(t, c, log_sigma, observed);
        } 
        
        // Calculate a part of the marginal.
        total = 0;
        for(size_t i=0; i<num_particles; ++i) 
            total += exp( W(t, i) );
        
        resampled = particles;
       
        for( size_t i=0; i<num_particles; ++i )
            p[i] = exp( W(t, i) ) / total;
        
        gsl_ran_discrete_t *g;
        g = gsl_ran_discrete_preproc(num_particles, p);
        
        for( size_t i=0; i<num_particles; ++i )
        {
            resample_particle_index = gsl_ran_discrete(r, g);
            particles[i] = resampled[resample_particle_index];
            std::cout << "Swapping particle " << i << " for particle " << resample_particle_index << std::endl;   
            std::cout << p[i] <<" against "<< p[resample_particle_index] <<std::endl;
        }
        
        /* This is working
        auto path = particles[i]->path();
        
        for( size_t k=0; k<K; ++k )
            for( size_t l=0; l<L; ++l )
            {
                 std::cout<< "Path at k="<<k<<"l="<<l<<"m="<<0<<" - "<<path(k,l,0,0)<<","<<path(k,l,0,1)<<std::endl;
                if(l<L-1)
                for( size_t m=1; m<M; ++m)
                {
                    std::cout<< "Path at k="<<k<<"l="<<l<<"m="<<m<<" - "<<path(k,l,m,0)<<","<<path(k,l,m,1)<<std::endl;
                }
               
            }
        */
    }
    
    std::string ts_file_name = "output/ts_smc.txt";
    std::ofstream ts_file;
    ts_file.open(ts_file_name);
    
    LangevinDynamics::PathType path;

    for(size_t l=0; l<L; ++l)
    {
        for(size_t k=0; k<K; ++k)
            for(size_t i=0; i<num_particles; i++)
            {
                path = particles[i]->path();
                ts_file << path(k,l,0,0) << "\t" << path(k,l,0,1) << "\t";
            }
            
        ts_file << std::endl;
        
        if(l<L-1)
        for(size_t m=0; m<M; ++m)
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
    
    ts_file_name = "output/ts_smc_observation.txt";
    ts_file.open(ts_file_name);
    
    for(size_t l=0; l<L; ++l)
    {
        for(size_t k=0; k<K; ++k)
            ts_file << observed(k,l,0) << "\t" << observed(k,l,1) << "\t";
        
        ts_file<<std::endl;
    }
    
    ts_file.close();
    
    particles.clear();
    resampled.clear();

}